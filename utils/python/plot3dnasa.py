#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python interface to NASA's PLOT3D file format with support for 3D
whole unformatted multi-block grid, solution and function files only.
Can handle both big-endian and little-endian byte ordering. Extends
NASA's original file format specification to include 64-bit record
lengths and quad-precision floating-point (as supported natively)."""

import sys
import numpy as np
from struct import pack, unpack_from
from itertools import izip

__all__ = ['Grid', 'Solution', 'Function', 'FileFormatError', 'fromfile',
           'cartesian_grid', 'cubic_bspline_support', 'tanh_support',
           'find_extents']

def fdcoeff(stencil, order=1):
    from scipy.linalg import solve
    from scipy.special import gamma
    p = np.arange(len(stencil))
    A = np.empty([p.size, p.size])
    for i, s in enumerate(stencil):
        A[i,:] = s ** p / gamma(p + 1)
    b = np.zeros(A.shape[1])
    b[order] = 1.
    return solve(A.T, b)

def fdmakeop(n, interior_accuracy=4, order=1,
             periodic=False, dtype=np.float64):
    from scipy.sparse import dok_matrix
    a = dok_matrix((n, n), dtype=dtype)
    m = (interior_accuracy + order - 1) / 2
    s = np.arange(-m, m + 1)
    if periodic:
        for k in xrange(n):
            c = fdcoeff(s, order)
            for ck, sk in izip(c, s):
                a[k,(k+sk)%n] = ck
        return a.tocsr()
    for k in xrange(m, n - m):
        c = fdcoeff(s, order)
        for ck, sk in izip(c, s):
            a[k,(k+sk)%n] = ck
    for k in range(m):
        s = np.arange(0, m + 1) - k
        c = fdcoeff(s, order)
        for ck, sk in izip(c, s):
            a[k,k+sk] = ck
            a[-(k+1),-(k+1+sk)] = (-1) ** order * ck
    return a.tocsr()

class FileFormatError(Exception):
    pass


class FileFormat(object):
    def __init__(self, filename=''):
        self.reclength_dtype = np.dtype(np.int32)
        self.integer_dtype = np.dtype(np.int32)
        self.real_dtype = np.dtype(np.float64)
        self.endianness = '='
        self.ncomponents = 0
        self.nblocks = 0
        self.size = np.empty([0, 3], dtype=np.dtype(np.int32))
        self.file_type = ''
        self.has_iblank = False
        self.aux_header = np.zeros([4], dtype=self.real_dtype)
        self.offsets = np.empty([0], dtype=np.dtype(np.int64))
        if filename:
            self.detect(filename)

    def detect(self, filename):
        try:
            with open(filename, 'rb') as f:
                # byte-ordering
                self._detect_endianness(f)
                if self.reclength_dtype is None:
                    raise FileFormatError(filename)
                # number of blocks
                f.seek(0)
                self.integer_dtype = self.integer_dtype.newbyteorder(
                    self.endianness)
                self.nblocks = self._get_integer(f)
                # block dimensions
                a = np.fromstring(self._get_raw_bytes(f), dtype=np.dtype(
                    np.int32).newbyteorder(self.endianness))
                if np.any(a <= 0):
                    raise FileFormatError(filename)
                if a.size == 4 * self.nblocks:
                    self.size = np.reshape(a, [4, self.nblocks],
                                           order = 'F')[:-1,:].T
                    self.ncomponents = a[3]
                else:
                    self.size = np.reshape(a, [3, self.nblocks], order = 'F').T
                # file type
                self._detect_file_type(f)
                if not self.file_type or self.real_dtype is None:
                    raise FileFormatError(filename)
                # aux header for solution files
                if self.file_type == 'solution':
                    self.aux_header = np.fromstring(self._get_raw_bytes(f),
                                                    self.real_dtype)
                # data offsets
                self._compute_offsets(f)
        except IOError:
            raise FileFormatError(filename)

    def _compute_offsets(self, f):
        self.offsets = np.empty(self.nblocks, dtype=np.dtype(np.int64))
        self.offsets[0] = f.tell() + self.reclength_dtype.itemsize
        for i in range(1, self.offsets.size):
            r = self._get_reclength(f)
            self.offsets[i] = self.offsets[i-1] + \
                              r + 2 * self.reclength_dtype.itemsize
            if i < self.offsets.size - 1 or self.file_type == 'solution':
                f.seek(r + 2 * self.reclength_dtype.itemsize, 1)
            if self.file_type == 'solution':
                r = self._get_reclength(f)
                self.offsets[i] += r + 2 * self.reclength_dtype.itemsize
                if i < self.offsets.size - 1:
                    f.seek(r + 2 * self.reclength_dtype.itemsize, 1)

    def _detect_endianness(self, f):
        s = f.read(8)
        for dtype in (endianness + reclength_dtype
                      for endianness in '=><' for reclength_dtype in 'iq'):
            if unpack_from(dtype, s)[0] == 4:
                self.reclength_dtype = np.dtype(dtype)
                self.endianness = self.reclength_dtype.str[0]
                return

    def _detect_file_type(self, f):
        r = self._get_reclength(f)
        if self.ncomponents > 0:
            self.file_type = 'function'
            r /= (self.ncomponents * np.prod(self.size[0,:]))
            self.real_dtype = np.dtype({4: np.float32, 8: np.float64,
                                        16: np.float128}[r])
            self.real_dtype = self.real_dtype.newbyteorder(self.endianness)
            return
        s = self._file_size(f) - f.tell()
        if r in [16, 32, 64]:
            if s == 2 * self.reclength_dtype.itemsize + r:
                self.file_type = 'grid'
                self.has_iblank = (np.prod(self.size[0,:]) == r / 16)
                self.real_dtype = np.dtype(np.float64)
                self.real_dtype = self.real_dtype.newbyteorder(self.endianness)
                return
            self.file_type = 'solution'
            self.real_dtype = np.dtype({16: np.float32, 32: np.float64,
                                        64: np.float128}[r])
            self.real_dtype = self.real_dtype.newbyteorder(self.endianness)
            return
        r /= np.prod(self.size[0,:])
        if r in [12, 16, 24, 28, 48, 52]:
            self.file_type = 'grid'
            self.has_iblank = r in [16, 28, 52]
            if self.has_iblank:
                r -= 4
            self.real_dtype = np.dtype({12: np.float32, 24: np.float64,
                                        48: np.float128}[r])
            self.real_dtype = self.real_dtype.newbyteorder(self.endianness)

    def _get_reclength(self, f, advance=False):
        n, = unpack_from(self.reclength_dtype.str,
                         f.read(self.reclength_dtype.itemsize))
        if not advance:
            f.seek(-self.reclength_dtype.itemsize, 1)
        return n

    def _get_raw_bytes(self, f):
        n, = unpack_from(self.reclength_dtype.str,
                         f.read(self.reclength_dtype.itemsize))
        s = f.read(n)
        f.seek(self.reclength_dtype.itemsize, 1)
        return s

    def _get_integer(self, f):
        return unpack_from(self.endianness + 'i',
                           f.read(2 * self.reclength_dtype.itemsize + 4),
                           offset = self.reclength_dtype.itemsize)[0]

    def _file_size(self, f):
        i = f.tell()
        f.seek(0, 2)
        s = f.tell()
        f.seek(i)
        return s


class MultiBlockCommon(object):
    def get_size(self, block_index=None):
        if block_index is None:
            return self.size
        return self.size[block_index,:]

    def set_size(self, size):
        ndim = np.array(size).ndim
        if ndim == 1:
            self.nblocks = 1
            self.size = np.empty([1, 3], dtype=np.dtype(
                    np.int32).newbyteorder(self._format.endianness))
            for i, s in enumerate(size):
                self.size[0,:] = [size[i] if i < len(size) else 1
                                  for i in range(3)]
        elif ndim == 2:
            self.nblocks = len(size)
            self.size = np.empty([self.nblocks, 3], dtype=np.dtype(
                    np.int32).newbyteorder(self._format.endianness))
            for i, s in enumerate(size):
                self.size[i,:] = [size[i][j] if j < len(size[i]) else 1
                                  for j in range(3)]

    def set_subzone(self, block_index, starts = None, ends = None):
        self._block_index = block_index
        self._subzone_starts, self._subzone_ends = None, None
        if block_index is None:
            self.set_size(self._format.size)
            return
        if block_index < 0 or block_index >= self._format.nblocks:
            raise ValueError('invalid block_index')
        self._subzone_starts = None
        self._subzone_ends = None
        size = self._format.size[block_index,:]
        if starts is not None and np.any(starts != 0):
            self._subzone_starts = np.array(starts)
        if ends is not None and (np.any(ends != -1) or
                                 np.any(ends != size - 1)):
            self._subzone_ends = np.array(ends)
            mask = self._subzone_ends < 0
            self._subzone_ends[mask] += size[mask]
        if self._subzone_starts is not None and self._subzone_ends is not None:
            if np.any(self._subzone_starts < 0) or \
               np.any(self._subzone_ends >= size) or \
               np.any(self._subzone_starts > self._subzone_ends):
                raise ValueError('invalid sub-zone')
            self.set_size([self._subzone_ends - self._subzone_starts + 1])
        else:
            self.set_size(self._format.size[block_index:block_index+1,:])
        return self

    def read_scalar(self, f, size, dtype, starts, ends):
        if starts is None or ends is None:
            return np.reshape(np.fromstring(
                f.read(dtype.itemsize * np.prod(size)), dtype),
                              size.tolist(), order='F')
        size_ = ends - starts + 1
        a = np.empty(size_.tolist(), dtype=dtype, order='F')
        f.seek(size[0] * size[1] * starts[2] * dtype.itemsize, 1)
        for k in range(starts[2], ends[2] + 1):
            f.seek(size[0] * starts[1] * dtype.itemsize, 1)
            for j in range(starts[1], ends[1] + 1):
                f.seek(starts[0] * dtype.itemsize, 1)
                a[:, j - starts[1], k - starts[2]] = np.fromstring(
                    f.read(size_[0] * dtype.itemsize), dtype=dtype)
                f.seek((size[0] - (ends[0] + 1)) * dtype.itemsize, 1)
            f.seek(size[0] * (size[1] - (ends[1] + 1)) * dtype.itemsize, 1)
        f.seek(size[0] * size[1] * (size[2] - (ends[2] + 1)) *
               dtype.itemsize, 1)
        return a

    def read_vector(self, f, size, dtype, n, starts, ends):
        if starts is None or ends is None:
            return np.reshape(np.fromstring(
                f.read(dtype.itemsize * n * np.prod(size)), dtype),
                              size.tolist() + [n], order='F')
        size_ = ends - starts + 1
        a = np.empty(size_.tolist() + [n], dtype=dtype, order='F')
        for i in range(n):
            f.seek(size[0] * size[1] * starts[2] * dtype.itemsize, 1)
            for k in range(starts[2], ends[2] + 1):
                f.seek(size[0] * starts[1] * dtype.itemsize, 1)
                for j in range(starts[1], ends[1] + 1):
                    f.seek(starts[0] * dtype.itemsize, 1)
                    a[:, j - starts[1], k - starts[2], i] = np.fromstring(
                          f.read(size_[0] * dtype.itemsize), dtype=dtype)
                    f.seek((size[0] - (ends[0] + 1)) * dtype.itemsize, 1)
                f.seek(size[0] * (size[1] - (ends[1] + 1)) * dtype.itemsize, 1)
            f.seek(size[0] * size[1] * (size[2] - (ends[2] + 1)) *
                   dtype.itemsize, 1)
        return a

    def write_header(self, f, ncomponents = 0):
        f.write(pack(self._format.reclength_dtype.str, 4))
        f.write(pack(self._format.endianness + 'i', self.nblocks))
        f.write(pack(self._format.reclength_dtype.str, 4))
        if ncomponents == 0:
            f.write(pack(self._format.reclength_dtype.str, 4 * self.size.size))
            f.write(self.size.tostring())
            f.write(pack(self._format.reclength_dtype.str, 4 * self.size.size))
        else:
            s = np.empty([self.nblocks, 4], dtype=np.dtype(
                    np.int32).newbyteorder(self._format.endianness))
            s[:,:-1] = self.size
            s[:,-1] = ncomponents
            f.write(pack(self._format.reclength_dtype.str, 4 * s.size))
            f.write(s.tostring())
            f.write(pack(self._format.reclength_dtype.str, 4 * s.size))

    def copy_from(self, a):
        self.filename = ''
        self._format = a._format
        return self.subzone_from(a)

    def subzone_from(self, a):
        self.set_subzone(a._block_index, a._subzone_starts, a._subzone_ends)
        self.set_size(a.get_size(), True)
        return self

class Grid(MultiBlockCommon):
    def __init__(self, filename='', block_index=None, subzone_starts=None,
                 subzone_ends=None, forceread=False):
        self.filename = filename
        self._format = FileFormat(filename)
        if self._format.file_type and self._format.file_type != 'grid':
            raise FileFormatError('%s is not a grid file' % filename)
        self.has_iblank = self._format.has_iblank
        self.set_subzone(block_index, subzone_starts, subzone_ends)
        if filename and forceread:
            self.load(filename)

    def set_size(self, size, allocate=False):
        super(Grid, self).set_size(size)
        self.xyz = [None] * self.nblocks
        if self.has_iblank:
            self.iblank = [None] * self.nblocks
        if allocate:
            self.allocate()
        return self

    def allocate(self):
        for i in range(self.nblocks):
            self.xyz[i] = np.empty(self.size[i].tolist() + [3],
                                   dtype=self._format.real_dtype, order='F')
        if self.has_iblank:
            self.iblank[i] = np.empty(self.size[i].tolist(),
                                      dtype=self._format.integer_dtype,
                                      order='F')
        return self

    def load(self, filename=''):
        if not filename:
            filename = self.filename
        elif not self.filename:
            self._format = FileFormat(filename)
            if self._format.file_type and self._format.file_type != 'grid':
                raise FileFormatError('%s is not a grid file' % filename)
            self.has_iblank = self._format.has_iblank
            self.set_subzone(None)
        with open(filename, 'rb') as f:
            j = 0
            for i in range(self._format.nblocks):
                if self._block_index is not None and i != self._block_index:
                    continue
                f.seek(self._format.offsets[i])
                self.xyz[j] = self.read_vector(
                    f, self._format.size[i,:], self._format.real_dtype, 3,
                    self._subzone_starts, self._subzone_ends)
                if self.has_iblank:
                    self.iblank[j] = self.read_scalar(
                        f, self._format.size[i,:], np.dtype(
                            np.int32).newbyteorder(self._format.endianness),
                        self._subzone_starts, self._subzone_ends)
                j += 1
        self.filename = filename
        return self

    def save(self, filename=''):
        if not filename:
            filename = self.filename
        with open(filename, 'wb') as f:
            self.write_header(f)
            for i in range(self.nblocks):
                s = 3 * self._format.real_dtype.itemsize
                if self.has_iblank:
                    s += self._format.integer_dtype.itemsize
                s *= np.prod(self.size[i,:])
                f.write(pack(self._format.reclength_dtype.str, s))
                f.write(self.xyz[i].tostring(order = 'F'))
                if self.has_iblank:
                    f.write(self.iblank[i].tostring(order = 'F'))
                f.write(pack(self._format.reclength_dtype.str, s))
        self.filename = filename

    def save_h5(self, filename=''):
        import h5py
        if not filename:
            filename = self.filename
        f = h5py.File(filename)
        d = dict(HEADER=np.zeros(1), numberOfGrids=self.nblocks)
        for key in d:
            f.attrs[key] = d[key]
        for i in range(self.nblocks):
            grp = f.require_group('Group%03d' % (i + 1))
            d = dict(HEADER=np.zeros(4), useIB=int(self.has_iblank),
                     gridSize=self.size[i,:], numberOfAuxVars=0)
            for key in d:
                grp.attrs[key] = d[key]
            for j in range(3):
                dset = grp.require_dataset('XYZ'[j], tuple(self.size[i,:]),
                                          dtype=self._format.real_dtype.str)
                dset[...] = self.xyz[i][:,:,:,j].flatten(order='F').reshape(
                    self.size[i,:], order='C')
            if self.has_iblank:
                dset = grp.require_dataset('IBLANK', tuple(self.size[i,:nd]),
                                          dtype=self._format.endianness + 'i')
                dset[...] = self.iblank[i].flatten(order='F').reshape(
                    self.size[i,:], order='C')
        f.close()

    def copy(self):
        g = Grid()
        g.filename = ''
        g._format = self._format
        g.has_iblank = self.has_iblank
        g.set_subzone(self._block_index, self._subzone_starts,
                      self._subzone_ends)
        g.set_size(self.get_size(), True)
        return g

    def minmax(self, block_index = None):
        if block_index is None:
            print 'Grid has %i block(s)\n' % self.nblocks
            for i in range(self.nblocks):
                print 'Block %i:' % (i + 1)
                self.minmax(i)
        else:
            labels = ['x', 'y', 'z']
            for i in range(3):
                print 'min. {0} = {1:+10.4E}, max. {0} = {2:+10.4E}'.format(
                    labels[i], self.xyz[block_index][:,:,:,i].min(),
                    self.xyz[block_index][:,:,:,i].max())
            print ''

    def __getitem__(self, key):
        return self.xyz[key]

    def __setitem__(self, key, item):
        self.xyz[key] = item

    def squeeze(self):
        for i, xyz in enumerate(self.xyz):
            nd = 1 if self.size[i,2] == 1 and self.size[i,1] == 1 else 2 \
                 if self.size[i,2] == 1 else 3
            self.xyz[i] = xyz[:,:,:,:nd]
        return self

class Solution(MultiBlockCommon):
    def __init__(self, filename='', block_index=None, subzone_starts=None,
                 subzone_ends=None, forceread=False):
        self.filename = filename
        self._format = FileFormat(filename)
        if self._format.aux_header is not None:
            self.time = self._format.aux_header[-1]
        else:
            self.time = 0.
        if self._format.file_type and self._format.file_type != 'solution':
            raise FileFormatError('%s is not a solution file' % filename)
        self.set_subzone(block_index, subzone_starts, subzone_ends)
        if filename and forceread:
            self.load(filename)

    def set_size(self, size, allocate=False):
        super(Solution, self).set_size(size)
        self.q = [None] * self.nblocks
        if allocate:
            self.allocate()
        return self

    def allocate(self):
        for i in range(self.nblocks):
            self.q[i] = np.empty(self.size[i].tolist() + [5],
                                 dtype=self._format.real_dtype, order='F')
        return self

    def load(self, filename=''):
        if not filename:
            filename = self.filename
        elif not self.filename:
            self._format = FileFormat(filename)
            if self._format.file_type and self._format.file_type != 'solution':
                raise FileFormatError('%s is not a solution file' % filename)
            if self._format.aux_header is not None:
                self.time = self._format.aux_header[-1]
            else:
                self.time = 0.
            self.set_subzone(None)
        with open(filename, 'rb') as f:
            f.seek(self._format.offsets[0] -
                   2 * self._format.reclength_dtype.itemsize -
                   4 * self._format.real_dtype.itemsize)
            self.time = np.fromstring(f.read(self._format.real_dtype.itemsize),
                                      dtype=self._format.real_dtype)[0]
            j = 0
            for i in range(self._format.nblocks):
                if self._block_index is not None and i != self._block_index:
                    continue
                f.seek(self._format.offsets[i])
                self.q[j] = self.read_vector(
                    f, self._format.size[i,:], self._format.real_dtype, 5,
                    self._subzone_starts, self._subzone_ends)
                j += 1
        self.filename = filename
        return self

    def save(self, filename=''):
        if not filename:
            filename = self.filename
        with open(filename, 'wb') as f:
            self.write_header(f)
            if self._format.aux_header is not None:
                aux_header = self._format.aux_header
            else:
                aux_header = np.zeros(dtype=self._format.real_dtype)
            aux_header[-1] = self.time
            for i in range(self.nblocks):
                s = 4 * self._format.real_dtype.itemsize
                f.write(pack(self._format.reclength_dtype.str, s))
                f.write(aux_header.tostring())
                f.write(pack(self._format.reclength_dtype.str, s))
                s = 5 * self._format.real_dtype.itemsize
                s *= np.prod(self.size[i,:])
                f.write(pack(self._format.reclength_dtype.str, s))
                f.write(self.q[i].tostring(order = 'F'))
                f.write(pack(self._format.reclength_dtype.str, s))
        self.filename = filename

    def save_h5(self, filename='', prefix='cv'):
        import h5py
        if not filename:
            filename = self.filename
        f = h5py.File(filename)
        d = dict(HEADER=np.zeros(1), numberOfGrids=self.nblocks)
        for key in d:
            f.attrs[key] = d[key]
        for i in range(self.nblocks):
            grp = f.require_group('Group%03d' % (i + 1))
            d = dict(HEADER=self._format.aux_header, gridSize=self.size[i,:],
                     numberOfAuxVars=0)
            for key in d:
                grp.attrs[key] = d[key]
            nd = 1 if self.size[i,2] == 1 and self.size[i,1] == 1 else 2 \
                 if self.size[i,2] == 1 else 3
            for k, j in enumerate([0] + range(1, nd + 1) + [4]):
                dset = grp.require_dataset('%s%02d' % (prefix, k + 1),
                                           tuple(self.size[i,:]),
                                           dtype=self._format.real_dtype.str)
                dset[...] = self.q[i][:,:,:,j].ravel(order='F').reshape(
                    self.size[i,:], order='C')
        f.close()

    def toprimitive(self, gamma = 1.4):
        for q in self.q:
            for i in range(1, 4):
                q[:,:,:,i] /= q[:,:,:,0]
            q[:,:,:,4] = (gamma - 1.) * (q[:,:,:,4] - 0.5 * q[:,:,:,0] *
                                         np.sum(q[:,:,:,1:4] ** 2, -1))
        return self

    def fromprimitive(self, gamma = 1.4):
        for q in self.q:
            q[:,:,:,4] /= gamma - 1.
            q[:,:,:,4] += 0.5 * q[:,:,:,0] * np.sum(q[:,:,:,1:4] ** 2, -1)
            for i in range(1, 4):
                q[:,:,:,i] *= q[:,:,:,0]
        return self

    def quiescent(self, gamma = 1.4):
        for q in self.q:
            q[:,:,:,0] = 1.
            q[:,:,:,1:4] = 0.
            q[:,:,:,4] = 1. / gamma
        return self

    def copy(self):
        s = Solution()
        s.filename = ''
        s._format = self._format
        s.time = self.time
        s.set_subzone(self._block_index, self._subzone_starts,
                      self._subzone_ends)
        s.set_size(self.get_size(), True)
        return s

    def minmax(self, block_index = None):
        if block_index is None:
            print 'Solution has %i block(s)\n' % self.nblocks
            for i in range(self.nblocks):
                print 'Block %i:' % (i + 1)
                self.minmax(i)
        else:
            labels = [u'ρ', u'ρu', u'ρv', u'ρw', u'e']
            for i in range(5):
                print u'min. {0:<2} = {1:+10.4E}, ' \
                    u'max. {0:<2} = {2:+10.4E}'.format(
                    labels[i], self.q[block_index][:,:,:,i].min(),
                    self.q[block_index][:,:,:,i].max())
            print ''

    def __getitem__(self, key):
        return self.q[key]

    def __setitem__(self, key, item):
        self.q[key] = item

    def squeeze(self):
        for i, q in enumerate(self.q):
            nd = 1 if self.size[i,2] == 1 and self.size[i,1] == 1 else 2 \
                 if self.size[i,2] == 1 else 3
            p = np.copy(q)
            self.q[i] = np.empty(self.size[i].tolist() + [nd + 2], order='F')
            for j in range(nd+1):
                self.q[i][:,:,:,j] = p[:,:,:,j]
            self.q[i][:,:,:,nd+1] = p[:,:,:,4]
        return self


class Function(MultiBlockCommon):
    def __init__(self, filename='', block_index=None, subzone_starts=None,
                 subzone_ends=None, forceread=False, ncomponents=1):
        self.filename = filename
        self._format = FileFormat(filename)
        if self._format.file_type and self._format.file_type != 'function':
            raise FileFormatError('%s is not a function file' % filename)
        if filename:
            self.ncomponents = self._format.ncomponents
        else:
            self.ncomponents = ncomponents
        self.set_subzone(block_index, subzone_starts, subzone_ends)
        if filename and forceread:
            self.load(filename)

    def set_size(self, size, allocate=False):
        super(Function, self).set_size(size)
        self.f = [None] * self.nblocks
        if allocate:
            self.allocate()
        return self

    def allocate(self):
        for i in range(self.nblocks):
            self.f[i] = np.empty(self.size[i].tolist() + [self.ncomponents],
                                 dtype=self._format.real_dtype, order='F')
        return self

    def load(self, filename=''):
        if not filename:
            filename = self.filename
        elif not self.filename:
            self._format = FileFormat(filename)
            if self._format.file_type and self._format.file_type != 'function':
                raise FileFormatError('%s is not a function file' % filename)
            self.ncomponents = self._format.ncomponents
            self.set_subzone(None)
        with open(filename, 'rb') as f:
            j = 0
            for i in range(self._format.nblocks):
                if self._block_index is not None and i != self._block_index:
                    continue
                f.seek(self._format.offsets[i])
                self.f[j] = self.read_vector(
                    f, self._format.size[i,:], self._format.real_dtype,
                    self.ncomponents, self._subzone_starts, self._subzone_ends)
                j += 1
        self.filename = filename
        return self

    def save(self, filename=''):
        if not filename:
            filename = self.filename
        with open(filename, 'wb') as f:
            self.write_header(f, self.ncomponents)
            for i in range(self.nblocks):
                s = self.ncomponents * self._format.real_dtype.itemsize
                s *= np.prod(self.size[i,:])
                f.write(pack(self._format.reclength_dtype.str, s))
                f.write(self.f[i].tostring(order = 'F'))
                f.write(pack(self._format.reclength_dtype.str, s))
        self.filename = filename

    def save_h5(self, filename='', prefix='f'):
        import h5py
        if not filename:
            filename = self.filename
        f = h5py.File(filename)
        d = dict(HEADER=np.zeros(1), numberOfGrids=self.nblocks)
        for key in d:
            f.attrs[key] = d[key]
        for i in range(self.nblocks):
            grp = f.require_group('Group%03d' % (i + 1))
            d = dict(HEADER=np.zeros(4), gridSize=self.size[i,:],
                     numberOfAuxVars=0)
            for key in d:
                grp.attrs[key] = d[key]
            nd = 1 if self.size[i,2] == 1 and self.size[i,1] == 1 else 2 \
                 if self.size[i,2] == 1 else 3
            for j in range(self.ncomponents):
                dset = grp.require_dataset('%s%02d' % (prefix, j + 1),
                                           tuple(self.size[i,:nd]),
                                           dtype=self._format.real_dtype.str)
                if nd == 3:
                    dset[...] = self.f[i][:,0,0,j]
                elif nd == 2:
                    dset[...] = self.f[i][:,:,0,j]
                else:
                    dset[...] = self.f[i][:,:,:,j]
        f.close()

    def copy(self):
        f = Function()
        f.filename = ''
        f._format = self._format
        f.ncomponents = self.ncomponents
        f.set_subzone(self._block_index, self._subzone_starts,
                      self._subzone_ends)
        f.set_size(self.get_size(), True)
        return f

    def __getitem__(self, key):
        return self.f[key]

    def __setitem__(self, key, item):
        self.f[key] = item


def fromfile(filename, block_index=None, subzone_starts=None,
             subzone_ends=None, file_type=''):
    file_type_ = file_type
    if not file_type_:
        file_type_ = FileFormat(filename).file_type
    if file_type_ == 'grid':
        return Grid(filename, block_index, subzone_starts, subzone_ends,
                    forceread = True)
    elif file_type_ == 'solution':
        return Solution(filename, block_index, subzone_starts, subzone_ends,
                        forceread = True)
    elif file_type_ == 'function':
        return Function(filename, block_index, subzone_starts, subzone_ends,
                        forceread = True)
    return None

def cartesian_grid(filename, block_index=0):
    f = FileFormat(filename)
    assert f.file_type == 'grid'
    g = Grid(filename, block_index)
    x = g.set_subzone(block_index, [0, 0, 0],
                      [-1, 0, 0]).load().xyz[0][:,0,0,0]
    y = g.set_subzone(block_index, [0, 0, 0],
                      [0, -1, 0]).load().xyz[0][0,:,0,1]
    z = g.set_subzone(block_index, [0, 0, 0],
                      [0, 0, -1]).load().xyz[0][0,0,:,2]
    return x, y, z

def cubic_bspline_support(x, x_min, x_max):
    from scipy.signal import cubic
    imin = np.argmin(np.abs(x - x_min))
    imax = np.argmin(np.abs(x - x_max))
    assert imax - imin + 1 >= 3
    return cubic(4. * (x - x[imin]) / (x[imax] - x[imin]) - 2.)

def tanh_support(x, x_min, x_max, sigma, xi):
    imin = np.argmin(np.abs(x - x_min))
    imax = np.argmin(np.abs(x - x_max))
    assert imax > imin
    f = lambda x: \
        np.tanh(sigma * (x + 1. - 0.5 * xi)) - \
        np.tanh(sigma * (x - 1. + 0.5 * xi))
    y = f(2. * (x - x[imin]) / (x[imax] - x[imin]) - 1.)
    return y - y.min()

def find_extents(x, x_min, x_max):
    return np.where(x >= x_min)[0][0], np.where(x <= x_max)[0][-1] + 2

def sbp(n, interior_accuracy=4, order=1, periodic=False, dtype=np.float64):
    if periodic:
        return fdmakeop(n, interior_accuracy, order, True, dtype)
    from scipy.sparse import dok_matrix
    assert interior_accuracy in [2, 4, 6, 8]
    assert order in [1, 2]
    a = dok_matrix((n, n), dtype=dtype)
    d = interior_accuracy if interior_accuracy != 2 else 1
    m = (interior_accuracy + order - 1) / 2
    s = np.arange(-m, m + 1)
    for k in xrange(d, n - d):
        c = fdcoeff(s, order)
        for ck, sk in izip(c, s):
            a[k,(k+sk)%n] = ck
    if order == 1:
        if interior_accuracy == 4:
            a[0,0:4] = a[-1,-1:-5:-1] = [-24./17., 59./34., -4./17., -3./34.]
            a[1,0:3] = a[-2,-1:-4:-1] = [-0.5, 0., 0.5]
            a[2,0:5] = a[-3,-1:-6:-1] = [4./43., -59./86., 0., 59./86.,
                                         -4./43.]
            a[3,0:6] = a[-4,-1:-7:-1] = [3./98., 0., -59./98., 0., 32./49,
                                         -4./49.]
    if order % 2 == 1:
        a[-d:,:] = -a[-d:,:]
    return a.tocsr()

def compute_jacobian(g, op_func=sbp, *op_func_args):
    f = Function(ncomponents=9).copy_from(g)
    for i, xyz in enumerate(g.xyz):
        n = g.get_size(i)
        f.f[i].fill(0.)
        for j in range(3):
            if n[j] > 1:
                op = op_func(n[j], *op_func_args)
                for k in range(3):
                    f.f[i][:,:,:,k+3*j] = np.apply_along_axis(
                        op.dot, j, xyz[:,:,:,k])
            else:
                f.f[i][:,:,:,j+3*j] = 1.
    return f

if __name__ == '__main__':
    pass
