"""Python interface to NASA's PLOT3D file format with support for 3D 
whole unformatted multi-block grid, solution and function files only. 
Can handle both big-endian and little-endian byte ordering. Extends 
NASA's original file format specification to include 64-bit record 
lengths and quad-precision floating-point (as supported natively)."""

import sys
import numpy as np
from struct import pack, unpack_from

__all__ = ['Grid', 'Solution', 'Function', 'FileFormatError', 'fromfile']

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
        self.size = np.empty([0, 3], dtype = np.dtype(np.int32))
        self.file_type = ''
        self.has_iblank = False
        self.aux_header = np.zeros([4], dtype = self.real_dtype)
        self.offsets = np.empty([0], dtype = np.dtype(np.int64))
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
                a = np.fromstring(self._get_raw_bytes(f),
                                  dtype = np.dtype(np.int32))
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
        self.offsets = np.empty(self.nblocks, dtype = np.dtype(np.int64))
        self.offsets[0] = f.tell() + self.reclength_dtype.itemsize
        for i in range(1, self.offsets.size):
            r = self._get_reclength(f)
            self.offsets[i] = self.offsets[i-1] + \
                              r + 2 * self.reclength_dtype.itemsize
            if i < self.offsets.size - 1:
                f.seek(r + 2 * self.reclength_dtype.itemsize, 1)
                if self.file_type == 'solution':
                    r = self._get_reclength(f)
                    self.offsets[i] += r + 2 * self.reclength_dtype.itemsize
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
            self.size = np.empty([1, 3], dtype = np.dtype(np.int32))
            for i, s in enumerate(size):
                self.size[0,:] = [size[i] if i < len(size) else 1
                                  for i in range(3)]            
        elif ndim == 2:
            self.nblocks = len(size)
            self.size = np.empty([self.nblocks, 3], dtype = np.dtype(np.int32))
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
        if ends is not None and np.any(ends != -1 and ends != size - 1):
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

    def read_scalar(self, f, size, dtype, starts, ends):
        if starts is None or ends is None:
            return np.reshape(np.fromstring(
                f.read(dtype.itemsize * np.prod(size)), dtype),
                              size.tolist(), order = 'F')
        size_ = ends - starts + 1
        a = np.empty(size_.tolist(), dtype = dtype)
        f.seek(size[0] * size[1] * starts[0] * dtype.itemsize, 1)
        for k in range(starts[2], ends[2] + 1):
            f.seek(size[0] * starts[0] * dtype.itemsize, 1)
            for j in range(starts[1], ends[1] + 1):
                f.seek(starts[0] * dtype.itemsize, 1)
                a[:, j - starts[1], k - starts[2]] = np.fromstring(
                    f.read(size_[0] * dtype.itemsize), dtype = dtype)
                f.seek((size[0] - (ends[0] + 1)) * dtype.itemsize, 1)
            f.seek(size[0] * (size[1] - (ends[1] + 1)) * dtype.itemsize, 1)
        f.seek(size[0] * size[1] * (size[2] - (ends[2] + 1)) *
               dtype.itemsize, 1)
        return a

    def read_vector(self, f, size, dtype, n, starts, ends):
        if starts is None or ends is None:
            return np.reshape(np.fromstring(
                f.read(dtype.itemsize * n * np.prod(size)), dtype),
                              size.tolist() + [n], order = 'F')
        size_ = ends - starts + 1
        a = np.empty(size_.tolist() + [n], dtype = dtype)
        for i in range(n):
            f.seek(size[0] * size[1] * starts[0] * dtype.itemsize, 1)
            for k in range(starts[2], ends[2] + 1):
                f.seek(size[0] * starts[0] * dtype.itemsize, 1)
                for j in range(starts[1], ends[1] + 1):
                    f.seek(starts[0] * dtype.itemsize, 1)
                    a[:, j - starts[1], k - starts[2], i] = np.fromstring(
                          f.read(size_[0] * dtype.itemsize), dtype = dtype)
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
            s = np.empty([self.nblocks, 4], dtype = np.dtype(np.int32))
            s[:,:-1] = self.size
            s[:,-1] = ncomponents
            f.write(pack(self._format.reclength_dtype.str, 4 * s.size))
            f.write(s.tostring())
            f.write(pack(self._format.reclength_dtype.str, 4 * s.size))

    def copy(self, a):
        self.filename = ''
        self._format = a._format
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
        if allocate:
            for i in range(self.nblocks):
                self.xyz[i] = np.empty(self.size[i].tolist() + [3], dtype =
                                       self._format.real_dtype)
        if self.has_iblank:
            self.iblank = [None] * self.nblocks
            if allocate:
                for i in range(self.nblocks):
                    self.iblank[i] = np.empty(self.size[i].tolist(), dtype =
                                              self._format.integer_dtype)
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
                        f, self._format.size[i,:], np.dtype(np.int32),
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

        
class Solution(MultiBlockCommon):
    def __init__(self, filename='', block_index=None, subzone_starts=None,
                 subzone_ends=None, forceread=False):
        self.filename = filename
        self._format = FileFormat(filename)
        if self._format.aux_header is not None:
            self.time = self._format.aux_header[0]
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
            for i in range(self.nblocks):
                self.q[i] = np.empty(self.size[i].tolist() + [5],
                                     dtype = self._format.real_dtype)
        return self

    def load(self, filename=''):
        if not filename:
            filename = self.filename
        elif not self.filename:
            self._format = FileFormat(filename)
            if self._format.file_type and self._format.file_type != 'solution':
                raise FileFormatError('%s is not a solution file' % filename)
            if self._format.aux_header is not None:
                self.time = self._format.aux_header[0]
            else:
                self.time = 0.            
            self.set_subzone(None)
        with open(filename, 'rb') as f:
            f.seek(self._format.offsets[0] -
                   2 * self._format.reclength_dtype.itemsize -
                   4 * self._format.real_dtype.itemsize)
            self.time = np.fromstring(f.read(self._format.real_dtype.itemsize),
                                      dtype = self._format.real_dtype)[0]
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
                aux_header = np.zeros(dtype = self._format.real_dtype)
            aux_header[0] = self.time
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
            for i in range(self.nblocks):
                self.f[i] = np.empty(self.size[i].tolist() +
                                     [self.ncomponents], dtype =
                                     self._format.real_dtype)
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

if __name__ == '__main__':
    pass