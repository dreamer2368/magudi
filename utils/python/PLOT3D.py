#!/usr/bin/env python
import numpy
import struct
try:
    from progressbar import ProgressBar
except:
    pass

class MultiBlockObject(object):

    def __init__(self):
        self._endianness = '='
        self._scalarType = 'd'
        self._offsetType = 'i'
        self._nGrids = 0
        self._size = []
        self._offsetBytes = 4
        self._scalarBytes = 8

    @property
    def endianness(self):
        return self._endianness

    @endianness.setter
    def endianness(self, value):
        allowed_values = ['=', '<', '>']
        if value not in allowed_values:
            raise ValueError("'%s' is not a valid endianness indicator (allowed values are: %s)" % (value, ", ".join(map(str, allowed_values))))
        self._endianness = value

    @property
    def scalarType(self):
        return self._scalarType

    @scalarType.setter
    def scalarType(self, value):
        allowed_values = ['d', 'f']
        if value not in allowed_values:
            raise ValueError("'%s' is not a valid format string for a scalar type (allowed values are: %s)" % (value, ", ".join(map(str, allowed_values))))
        self._scalarType = value
        if value == 'd':
            self._scalarBytes = 8
        if value == 'f':
            self._scalarBytes = 4

    @property
    def offsetType(self):
        return self._offsetType

    @offsetType.setter
    def offsetType(self, value):
        allowed_values = ['i', 'q']
        if value not in allowed_values:
            raise ValueError("'%s' is not a valid format string for an offset type (allowed values are: %s)" % (value, ", ".join(map(str, allowed_values))))
        self._offsetType = value
        if value == 'i':
            self._offsetBytes = 4
        if value == 'q':
            self._offsetBytes = 8

    @property
    def scalarBytes(self):
        return self._scalarBytes

    @property
    def offsetBytes(self):
        return self._offsetBytes                

    @property
    def nGrids(self):
        return self._nGrids

    @nGrids.setter
    def nGrids(self, value):
        if not isinstance(value, (int, long)):
            raise TypeError("Number of grids must be an integer")        
        if value <= 0:
            raise ValueError("Number of grids must be positive")
        self._nGrids = value
        self._size = [ [0, 0, 0] for iGrid in range(0, self._nGrids) ]
        self.Update()

    @property
    def size(self):
        return self._size

    def GetSize(self, grid_index = None):
        if grid_index is None:
            return self._size
        if not isinstance(grid_index, (int, long)):
            raise TypeError("Grid index must be an integer")
        if grid_index < 0 or grid_index >= self._nGrids:
            raise IndexError("Grid index out of range (expected a number in [%i,%i])" % (0, self._nGrids - 1))
        return self._size[grid_index]

    def SetSize(self, grid_index, size, allocate = True):
        if not isinstance(grid_index, (int, long)):
            raise TypeError("Grid index must be an integer")
        if grid_index < 0 or grid_index >= self._nGrids:
            raise IndexError("Grid index out of range (expected a number in [%i,%i])" % (0, self._nGrids - 1))
        assert len(size) == 3
        for idim in range(0, 3):
            self._size[grid_index][idim] = size[idim]
        if allocate is True:
            self.Update(grid_index)
        return None

    def GetNumberOfPoints(self):
        number_of_points = 0
        for iGrid in range(0, self._nGrids):
            number_of_points = number_of_points + numpy.product(self._size[iGrid])
        return number_of_points

    def CopyFrom(self, obj):
        self.endianness = obj.endianness
        self.scalarType = obj.scalarType
        self.offsetType = obj.offsetType
        if obj.nGrids > 0:
            self.nGrids = obj.nGrids
            for iGrid in range(0, self._nGrids):
                self.SetSize(iGrid, obj.GetSize(iGrid))
                self.Update(iGrid)
        return None

class Grid(MultiBlockObject):

    def __init__(self):
        super(Grid, self).__init__()
        self._hasIBLANK = True
        self._X = []
        self._IBLANK = []

    @property
    def hasIBLANK(self):
        return self._hasIBLANK

    @hasIBLANK.setter
    def hasIBLANK(self, value):
        assert isinstance(value, bool)
        self._hasIBLANK = value

    @property
    def X(self):
        return self._X

    @property
    def IBLANK(self):
        return self._IBLANK

    def Update(self, grid_index = None):
        if grid_index is None:
            self._X = [ None for iGrid in range(0, self.nGrids) ]
            self._IBLANK = [ None for iGrid in range(0, self.nGrids) ]
        else:
            grid_size = self.GetSize(grid_index)
            self._X[grid_index] = numpy.zeros(grid_size + [3], dtype = self.endianness + self.scalarType)
            self._IBLANK[grid_index] = numpy.ones(grid_size, dtype = self.endianness + 'i')
        return None

    def Export(self, filename):
        f = open(filename, "wb")
        record_size = 4
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        f.write(struct.pack(self.endianness + 'i', self.nGrids))
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        record_size = 12 * self.nGrids
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        for iGrid in range(0, self.nGrids):
            f.write(struct.pack(self.endianness + '3i', *self.GetSize(iGrid)))
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        for iGrid in range(0, self.nGrids):
            record_size = (4 + 3 * self.scalarBytes) * numpy.product(self.GetSize(iGrid))
            f.write(struct.pack(self.endianness + self.offsetType, record_size))
            f.write(self._X[iGrid].tostring(order = 'F'))
            f.write(self._IBLANK[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.endianness + self.offsetType, record_size))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        self.nGrids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, struct.unpack(self.endianness + '3i', f.read(12)))
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(self.offsetBytes, 1)
            self._X[iGrid][:,:,:,:] = numpy.reshape(struct.unpack(self.endianness + str(3 * numpy.product(self.GetSize(iGrid))) + self.scalarType, f.read(3 * self.scalarBytes * numpy.product(self.GetSize(iGrid)))), self.GetSize(iGrid) + [3], order = 'F')
            if self._hasIBLANK is True:
                self._IBLANK[iGrid][:,:,:] = numpy.reshape(struct.unpack(self.endianness + str(numpy.product(self.GetSize(iGrid))) + 'i', f.read(4 * numpy.product(self.GetSize(iGrid)))), self.GetSize(iGrid), order = 'F')
            f.seek(self.offsetBytes, 1)
        f.close()
        return None

    def ImportSkeleton(self, filename):
        f = open(filename, "rb")
        self.nGrids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, struct.unpack(self.endianness + '3i', f.read(12)), allocate = False)
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(2 * self.offsetBytes + 3 * self.scalarBytes * numpy.product(self.GetSize(iGrid)), 1)
            if self._hasIBLANK is True:
                f.seek(4 * numpy.product(self.GetSize(iGrid)), 1)
        f.close()
        return None

    def ImportSubarray(self, filename, grid_index, start_indices, end_indices, show_progress = False):

        f = open(filename, "rb")
        number_of_grids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        grid_size = [ [0, 0, 0] for iGrid in range(0, number_of_grids) ]
        self.nGrids = 1
 
        assert len(start_indices) == 3 and len(end_indices) == 3       
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, number_of_grids):
            grid_size[iGrid] = struct.unpack(self.endianness + '3i', f.read(12))
            if iGrid == grid_index:
                for idim in range(0, 3):
                    if start_indices[idim] < 0:
                        start_indices[idim] = start_indices[idim] + grid_size[iGrid][idim]
                    if end_indices[idim] < 0:
                        end_indices[idim] = end_indices[idim] + grid_size[iGrid][idim]
                    assert start_indices[idim] >= 0 and end_indices[idim] < grid_size[iGrid][idim] and end_indices[idim] >= start_indices[idim]
        f.seek(self.offsetBytes, 1)
        self.SetSize(0, [ end_indices[idim] - start_indices[idim] + 1 for idim in range(0, 3) ])
        
        for iGrid in range(0, number_of_grids):
            f.seek(self.offsetBytes, 1)
            if iGrid != grid_index:
                f.seek(3 * self.scalarBytes * numpy.product(grid_size[iGrid]), 1)
                if self._hasIBLANK is True:
                    f.seek(4 * numpy.product(grid_size[iGrid]), 1)
            else:
                
                if show_progress is True:
                    pbar = ProgressBar(maxval = 3 * (end_indices[2] - start_indices[2] + 1) * (end_indices[1] - start_indices[1] + 1))
                    pbar.start()
                for idim in range(0, 3):
                    f.seek(self.scalarBytes * grid_size[iGrid][0] * grid_size[iGrid][1] * start_indices[2], 1)
                    for k in range(start_indices[2], end_indices[2] + 1):
                        f.seek(self.scalarBytes * grid_size[iGrid][0] * start_indices[1], 1)
                        for j in range(start_indices[1], end_indices[1] + 1):
                            if show_progress is True:
                                pbar.update(j - start_indices[1] + 1 + (end_indices[1] - start_indices[1] + 1) * (k - start_indices[2] + (end_indices[2] - start_indices[2] + 1) * idim))
                            f.seek(self.scalarBytes * start_indices[0], 1)
                            self._X[0][:, j - start_indices[1], k - start_indices[2], idim] = numpy.array(struct.unpack(self.endianness + str((end_indices[0] - start_indices[0] + 1)) + self.scalarType, f.read(self.scalarBytes * (end_indices[0] - start_indices[0] + 1))), dtype = self.endianness + self.scalarType)
                            f.seek(self.scalarBytes * (grid_size[iGrid][0] - (end_indices[0] + 1)), 1)
                        f.seek(self.scalarBytes * grid_size[iGrid][0] * (grid_size[iGrid][1] - (end_indices[1] + 1)), 1)
                    f.seek(self.scalarBytes * grid_size[iGrid][0] * grid_size[iGrid][1] * (grid_size[iGrid][2] - (end_indices[2] + 1)), 1)
                if show_progress is True:
                    pbar.finish()

                if self._hasIBLANK is True:
                    if show_progress is True:
                        pbar = ProgressBar(maxval = (end_indices[2] - start_indices[2] + 1) * (end_indices[1] - start_indices[1] + 1))
                        pbar.start()
                    f.seek(4 * grid_size[iGrid][0] * grid_size[iGrid][1] * start_indices[2], 1)
                    for k in range(start_indices[2], end_indices[2] + 1):
                        f.seek(4 * grid_size[iGrid][0] * start_indices[1], 1)
                        for j in range(start_indices[1], end_indices[1] + 1):
                            if show_progress is True:
                                pbar.update(j - start_indices[1] + 1 + (end_indices[1] - start_indices[1] + 1) * (k - start_indices[2]))
                            f.seek(4 * start_indices[0], 1)
                            self._IBLANK[0][:, j - start_indices[1], k - start_indices[2]] = numpy.array(struct.unpack(self.endianness + str((end_indices[0] - start_indices[0] + 1)) + 'i', f.read(4 * (end_indices[0] - start_indices[0] + 1))), dtype = self.endianness + 'i')
                            f.seek(4 * (grid_size[iGrid][0] - (end_indices[0] + 1)), 1)
                        f.seek(4 * grid_size[iGrid][0] * (grid_size[iGrid][1] - (end_indices[1] + 1)), 1)
                    f.seek(4 * grid_size[iGrid][0] * grid_size[iGrid][1] * (grid_size[iGrid][2] - (end_indices[2] + 1)), 1)
                    if show_progress is True:
                        pbar.finish()
                    
            f.seek(self.offsetBytes, 1)
        f.close()

        return None

class Solution(MultiBlockObject):

    def __init__(self):
        super(Solution, self).__init__()
        self._Q = []

    @property
    def Q(self):
        return self._Q

    def Update(self, grid_index = None):
        if grid_index is None:
            self._Q = [ None for iGrid in range(0, self.nGrids) ]
        else:
            grid_size = self.GetSize(grid_index)
            self._Q[grid_index] = numpy.zeros(grid_size + [5], dtype = self.endianness + self.scalarType)
        return None

    def Export(self, filename, Q_header = None):
        f = open(filename, "wb")
        record_size = 4
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        f.write(struct.pack(self.endianness + 'i', self.nGrids))
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        record_size = 12 * self.nGrids
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        for iGrid in range(0, self.nGrids):
            f.write(struct.pack(self.endianness + '3i', *self.GetSize(iGrid)))
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        if Q_header is None:
            Q_header = numpy.zeros(4, dtype = self.endianness + self.scalarType)        
        for iGrid in range(0, self.nGrids):
            record_size = 4 * self.scalarBytes
            f.write(struct.pack(self.endianness + self.offsetType, record_size))
            f.write(Q_header.tostring())
            f.write(struct.pack(self.endianness + self.offsetType, record_size))
            record_size = 5 * self.scalarBytes * numpy.product(self.GetSize(iGrid))
            f.write(struct.pack(self.endianness + self.offsetType, record_size))
            f.write(self._Q[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.endianness + self.offsetType, record_size))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        self.nGrids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, struct.unpack(self.endianness + '3i', f.read(12)))
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(3 * self.offsetBytes + 4 * self.scalarBytes, 1)
            self._Q[iGrid][:,:,:,:] = numpy.reshape(struct.unpack(self.endianness + str(5 * numpy.product(self.GetSize(iGrid))) + self.scalarType, f.read(5 * self.scalarBytes * numpy.product(self.GetSize(iGrid)))), self.GetSize(iGrid) + [5], order = 'F')
            f.seek(self.offsetBytes, 1)
        f.close()
        return None

    def ImportSkeleton(self, filename):
        f = open(filename, "rb")
        self.nGrids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, struct.unpack(self.endianness + '3i', f.read(12)), allocate = False)
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(4 * self.offsetBytes + (4 + 5 * numpy.product(self.GetSize(iGrid))) * self.scalarBytes, 1)
        f.close()
        return None

    def ImportSubarray(self, filename, grid_index, start_indices, end_indices, show_progress = False):

        f = open(filename, "rb")
        number_of_grids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        grid_size = [ [0, 0, 0] for iGrid in range(0, number_of_grids) ]
        self.nGrids = 1
 
        assert len(start_indices) == 3 and len(end_indices) == 3       
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, number_of_grids):
            grid_size[iGrid] = struct.unpack(self.endianness + '3i', f.read(12))
            if iGrid == grid_index:
                for idim in range(0, 3):
                    if start_indices[idim] < 0:
                        start_indices[idim] = start_indices[idim] + grid_size[iGrid][idim]
                    if end_indices[idim] < 0:
                        end_indices[idim] = end_indices[idim] + grid_size[iGrid][idim]
                    assert start_indices[idim] >= 0 and end_indices[idim] < grid_size[iGrid][idim] and end_indices[idim] >= start_indices[idim]
        f.seek(self.offsetBytes, 1)
        self.SetSize(0, [ end_indices[idim] - start_indices[idim] + 1 for idim in range(0, 3) ])
        
        for iGrid in range(0, number_of_grids):
            f.seek(3 * self.offsetBytes + 4 * self.scalarBytes, 1)
            if iGrid != grid_index:
                f.seek(5 * numpy.product(grid_size[iGrid]) * self.scalarBytes, 1)
            else:
                if show_progress is True:
                    pbar = ProgressBar(maxval = 5 * (end_indices[2] - start_indices[2] + 1) * (end_indices[1] - start_indices[1] + 1))
                    pbar.start()
                for l in range(0, 5):
                    f.seek(self.scalarBytes * grid_size[iGrid][0] * grid_size[iGrid][1] * start_indices[2], 1)
                    for k in range(start_indices[2], end_indices[2] + 1):
                        f.seek(self.scalarBytes * grid_size[iGrid][0] * start_indices[1], 1)
                        for j in range(start_indices[1], end_indices[1] + 1):
                            if show_progress is True:
                                pbar.update(j - start_indices[1] + 1 + (end_indices[1] - start_indices[1] + 1) * (k - start_indices[2] + (end_indices[2] - start_indices[2] + 1) * l))                            
                            f.seek(self.scalarBytes * start_indices[0], 1)
                            self._Q[0][:, j - start_indices[1], k - start_indices[2], l] = numpy.array(struct.unpack(self.endianness + str((end_indices[0] - start_indices[0] + 1)) + self.scalarType, f.read(self.scalarBytes * (end_indices[0] - start_indices[0] + 1))), dtype = self.endianness + self.scalarType)
                            f.seek(self.scalarBytes * (grid_size[iGrid][0] - (end_indices[0] + 1)), 1)
                        f.seek(self.scalarBytes * grid_size[iGrid][0] * (grid_size[iGrid][1] - (end_indices[1] + 1)), 1)
                    f.seek(self.scalarBytes * grid_size[iGrid][0] * grid_size[iGrid][1] * (grid_size[iGrid][2] - (end_indices[2] + 1)), 1)
                if show_progress is True:
                    pbar.finish()
            f.seek(self.offsetBytes, 1)
        f.close()

        return None

class Function(MultiBlockObject):

    def __init__(self):
        super(Function, self).__init__()
        self._nComponents = 1
        self._F = []

    @property
    def nComponents(self):
        return self._nComponents

    @nComponents.setter
    def nComponents(self, value):
        if not isinstance(value, (int, long)):
            raise TypeError("Number of components must be an integer")
        if value < 0:
            raise IndexError("Number of components must be positive")
        self._nComponents = value

    @property
    def F(self):
        return self._F

    def Update(self, grid_index = None):
        if grid_index is None:
            self._F = [ None for iGrid in range(0, self.nGrids) ]
        else:
            grid_size = self.GetSize(grid_index)
            self._F[grid_index] = numpy.zeros(grid_size + [self._nComponents], dtype = self.endianness + self.scalarType)
        return None

    def Export(self, filename, Q_header = None):
        f = open(filename, "wb")
        record_size = 4
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        f.write(struct.pack(self.endianness + 'i', self.nGrids))
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        record_size = 16 * self.nGrids
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        for iGrid in range(0, self.nGrids):
            f.write(struct.pack(self.endianness + '3i', *self.GetSize(iGrid)))
            f.write(struct.pack(self.endianness + 'i', self._nComponents))
        f.write(struct.pack(self.endianness + self.offsetType, record_size))
        for iGrid in range(0, self.nGrids):
            record_size = self._nComponents * self.scalarBytes * numpy.product(self.GetSize(iGrid))
            f.write(struct.pack(self.endianness + self.offsetType, record_size))
            f.write(self._F[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.endianness + self.offsetType, record_size))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        self.nGrids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, struct.unpack(self.endianness + '3i', f.read(12)), allocate = False)
            if iGrid == 0:
                self._nComponents, = struct.unpack(self.endianness + 'i', f.read(4))
            else:
                number_of_components, = struct.unpack(self.endianness + 'i', f.read(4))
                if number_of_components != self._nComponents:
                    raise IOError("%s: Invalid PLOT3D function file: number of components in block %i (= %i) differs from %i" % (filename, iGrid + 1, number_of_components, self._nComponents))
            self.Update(iGrid)
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(self.offsetBytes, 1)
            self._F[iGrid][:,:,:,:] = numpy.reshape(struct.unpack(self.endianness + str(self._nComponents * numpy.product(self.GetSize(iGrid))) + self.scalarType, f.read(self._nComponents * self.scalarBytes * numpy.product(self.GetSize(iGrid)))), self.GetSize(iGrid) + [self._nComponents], order = 'F')
            f.seek(self.offsetBytes, 1)
        f.close()
        return None

    def ImportSkeleton(self, filename):
        f = open(filename, "rb")
        self.nGrids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, struct.unpack(self.endianness + '3i', f.read(12)), allocate = False)
            if iGrid == 0:
                self._nComponents, = struct.unpack(self.endianness + 'i', f.read(4))
            else:
                number_of_components, = struct.unpack(self.endianness + 'i', f.read(4))
                if number_of_components != self._nComponents:
                    raise IOError("%s: Invalid PLOT3D function file: number of components in block %i (= %i) differs from %i" % (filename, iGrid + 1, number_of_components, self._nComponents))
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(2 * self.offsetBytes + (self._nComponents * numpy.product(self.GetSize(iGrid))) * self.scalarType, 1)
        f.close()
        return None

    def ImportSubarray(self, filename, grid_index, start_indices, end_indices, show_progress = False):

        f = open(filename, "rb")
        number_of_grids, = struct.unpack(self.endianness + 'i', f.read(2 * self.offsetBytes + 4)[self.offsetBytes:-self.offsetBytes])
        grid_size = [ [0, 0, 0] for iGrid in range(0, number_of_grids) ]
        self.nGrids = 1
 
        assert len(start_indices) == 3 and len(end_indices) == 3       
        f.seek(self.offsetBytes, 1)
        for iGrid in range(0, number_of_grids):
            grid_size[iGrid] = struct.unpack(self.endianness + '3i', f.read(12))
            if iGrid == 0:
                self._nComponents, = struct.unpack(self.endianness + 'i', f.read(4))
            else:
                number_of_components, = struct.unpack(self.endianness + 'i', f.read(4))
                if number_of_components != self._nComponents:
                    raise IOError("%s: Invalid PLOT3D function file: number of components in block %i (= %i) differs from %i" % (filename, iGrid + 1, number_of_components, self._nComponents))
            if iGrid == grid_index:
                for idim in range(0, 3):
                    if start_indices[idim] < 0:
                        start_indices[idim] = start_indices[idim] + grid_size[iGrid][idim]
                    if end_indices[idim] < 0:
                        end_indices[idim] = end_indices[idim] + grid_size[iGrid][idim]
                    assert start_indices[idim] >= 0 and end_indices[idim] < grid_size[iGrid][idim] and end_indices[idim] >= start_indices[idim]
        f.seek(self.offsetBytes, 1)
        self.SetSize(0, [ end_indices[idim] - start_indices[idim] + 1 for idim in range(0, 3) ])
        
        for iGrid in range(0, number_of_grids):
            f.seek(self.offsetBytes, 1)
            if iGrid != grid_index:
                f.seek(self._nComponents * numpy.product(grid_size[iGrid]) * self.scalarBytes, 1)
            else:
                if show_progress is True:
                    pbar = ProgressBar(maxval = self._nComponents * (end_indices[2] - start_indices[2] + 1) * (end_indices[1] - start_indices[1] + 1))
                    pbar.start()                
                for l in range(0, self._nComponents):
                    f.seek(self.scalarBytes * grid_size[iGrid][0] * grid_size[iGrid][1] * start_indices[2], 1)
                    for k in range(start_indices[2], end_indices[2] + 1):
                        f.seek(self.scalarBytes * grid_size[iGrid][0] * start_indices[1], 1)
                        for j in range(start_indices[1], end_indices[1] + 1):
                            if show_progress is True:
                                pbar.update(j - start_indices[1] + 1 + (end_indices[1] - start_indices[1] + 1) * (k - start_indices[2] + (end_indices[2] - start_indices[2] + 1) * l))
                            f.seek(self.scalarBytes * start_indices[0], 1)
                            self._F[0][:, j - start_indices[1], k - start_indices[2], l] = numpy.array(struct.unpack(self.endianness + str((end_indices[0] - start_indices[0] + 1)) + self.scalarType, f.read(self.scalarBytes * (end_indices[0] - start_indices[0] + 1))), dtype = self.endianness + self.scalarType)
                            f.seek(self.scalarBytes * (grid_size[iGrid][0] - (end_indices[0] + 1)), 1)
                        f.seek(self.scalarBytes * grid_size[iGrid][0] * (grid_size[iGrid][1] - (end_indices[1] + 1)), 1)
                    f.seek(self.scalarBytes * grid_size[iGrid][0] * grid_size[iGrid][1] * (grid_size[iGrid][2] - (end_indices[2] + 1)), 1)
                if show_progress is True:
                    pbar.finish()
            f.seek(self.offsetBytes, 1)
        f.close()

        return None    
