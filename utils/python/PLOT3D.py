#!/usr/bin/env python
import numpy as np
import struct
try:
    from progressbar import ProgressBar
except:
    pass

class MultiBlockObject(object):

    def __init__(self):
        self._endianness = '='
        self._integerType = np.dtype(np.int32)
        self._scalarType = np.dtype(np.float64)
        self._offsetType = np.dtype(np.int32)
        self._nGrids = 0
        self._size = []
        self._integerTypeStr = '=i'
        self._offsetTypeStr = '=i'

    @property
    def endianness(self):
        return self._endianness

    @endianness.setter
    def endianness(self, value):
        allowed_values = ['=', '<', '>']
        if value not in allowed_values:
            raise ValueError("'%s' is not a valid endianness indicator (allowed values are: %s)" % (value, ", ".join(map(str, allowed_values))))
        self._endianness = value
        self._integerType = self._integerType.newbyteorder(self._endianness)
        self._scalarType = self._scalarType.newbyteorder(self._endianness)
        self._offsetType = self._offsetType.newbyteorder(self._endianness)
        self._integerTypeStr = self._endianness + self._integerTypeStr[1:]
        self._offsetTypeStr = self._endianness + self._offsetTypeStr[1:]

    @property
    def scalarType(self):
        return self._scalarType

    @scalarType.setter
    def scalarType(self, value):
        allowed_values = np.sctypes['float'] + np.sctypes['complex']
        if value not in allowed_values:
            raise ValueError("'%s' is not a valid format string for a scalar type (allowed values are: %s)" % (value, ", ".join(map(str, allowed_values))))
        self._scalarType = np.dtype(value).newbyteorder(self._endianness)

    @property
    def integerType(self):
        return self._integerType

    @integerType.setter
    def integerType(self, value):
        allowed_values = np.sctypes['int'] + np.sctypes['uint']
        if value not in allowed_values:
            raise ValueError("'%s' is not a valid format string for a integer type (allowed values are: %s)" % (value, ", ".join(map(str, allowed_values))))
        self._integerType = np.dtype(value).newbyteorder(self._endianness)
        if value == np.uint8:
            self._integerTypeStr = self._endianness + 'B'
        elif value == np.uint16:
            self._integerTypeStr = self._endianness + 'H'
        elif value == np.uint32:
            self._integerTypeStr = self._endianness + 'I'
        elif value == np.uint64:
            self._integerTypeStr = self._endianness + 'Q'
        elif value == np.int8:
            self._integerTypeStr = self._endianness + 'b'
        elif value == np.int16:
            self._integerTypeStr = self._endianness + 'h'
        elif value == np.int32:
            self._integerTypeStr = self._endianness + 'i'
        elif value == np.int64:
            self._integerTypeStr = self._endianness + 'q'

    @property
    def offsetType(self):
        return self._offsetType

    @offsetType.setter
    def offsetType(self, value):
        allowed_values = np.sctypes['uint']
        if value not in allowed_values:
            raise ValueError("'%s' is not a valid format string for an offset type (allowed values are: %s)" % (value, ", ".join(map(str, allowed_values))))
        self._offsetType = np.dtype(value).newbyteorder(self._endianness)
        if value == np.uint8:
            self._offsetTypeStr = self._endianness + 'B'
        elif value == np.uint16:
            self._offsetTypeStr = self._endianness + 'H'
        elif value == np.uint32:
            self._offsetTypeStr = self._endianness + 'I'
        elif value == np.uint64:
            self._offsetTypeStr = self._endianness + 'Q'
        elif value == np.int8:
            self._offsetTypeStr = self._endianness + 'b'
        elif value == np.int16:
            self._offsetTypeStr = self._endianness + 'h'
        elif value == np.int32:
            self._offsetTypeStr = self._endianness + 'i'
        elif value == np.int64:
            self._offsetTypeStr = self._endianness + 'q'

    @property
    def offsetTypeStr(self):
        return self._offsetTypeStr

    @property
    def integerTypeStr(self):
        return self._integerTypeStr

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
            number_of_points = number_of_points + np.product(self._size[iGrid])
        return number_of_points

    def CopyFrom(self, obj):
        self._endianness = obj.endianness
        self._integerType = obj.integerType
        self._scalarType = obj.scalarType
        self._offsetType = obj.offsetType
        self._integerTypeStr = obj.integerTypeStr
        self._offsetTypeStr = obj.offsetTypeStr
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
            self._X[grid_index] = np.zeros(grid_size + [3], dtype = self.scalarType)
            self._IBLANK[grid_index] = np.ones(grid_size, dtype = self.integerType)
        return None

    def Export(self, filename):
        f = open(filename, "wb")
        record_size = self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, record_size))
        f.write(struct.pack(self.integerTypeStr, self.nGrids))
        f.write(struct.pack(self.offsetTypeStr, record_size))
        record_size = 3 * self.nGrids * self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, record_size))
        for iGrid in range(0, self.nGrids):
            gridSize = np.array(self.GetSize(iGrid), dtype = self.integerType)
            f.write(gridSize.tostring())
        f.write(struct.pack(self.offsetTypeStr, record_size))
        for iGrid in range(0, self.nGrids):
            record_size = (self.integerType.itemsize + 3 * self.scalarType.itemsize) * np.product(self.GetSize(iGrid))
            f.write(struct.pack(self.offsetTypeStr, record_size))
            f.write(self._X[iGrid].tostring(order = 'F'))
            f.write(self._IBLANK[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.offsetTypeStr, record_size))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType))
        f.seek(self.offsetType.itemsize, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(self.offsetType.itemsize, 1)
            dtype = 3 * np.product(self.GetSize(iGrid)) * self.scalarType
            self._X[iGrid][:,:,:,:] = np.reshape(np.fromstring(f.read(dtype.itemsize), dtype = dtype), self.GetSize(iGrid) + [3], order = 'F')
            if self._hasIBLANK is True:
                dtype = np.product(self.GetSize(iGrid)) * self.integerType
                self._IBLANK[iGrid][:,:,:] = np.reshape(np.fromstring(f.read(dtype.itemsize), dtype = dtype), self.GetSize(iGrid), order = 'F')
            f.seek(self.offsetType.itemsize, 1)
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
            f.seek(2 * self.offsetBytes + 3 * self.scalarBytes * np.product(self.GetSize(iGrid)), 1)
            if self._hasIBLANK is True:
                f.seek(4 * np.product(self.GetSize(iGrid)), 1)
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
                f.seek(3 * self.scalarBytes * np.product(grid_size[iGrid]), 1)
                if self._hasIBLANK is True:
                    f.seek(4 * np.product(grid_size[iGrid]), 1)
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
                            self._X[0][:, j - start_indices[1], k - start_indices[2], idim] = np.array(struct.unpack(self.endianness + str((end_indices[0] - start_indices[0] + 1)) + self.scalarType, f.read(self.scalarBytes * (end_indices[0] - start_indices[0] + 1))), dtype = self.scalarType.newbyteorder(self.endianness))
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
                            self._IBLANK[0][:, j - start_indices[1], k - start_indices[2]] = np.array(struct.unpack(self.endianness + str((end_indices[0] - start_indices[0] + 1)) + 'i', f.read(4 * (end_indices[0] - start_indices[0] + 1))), dtype = np.int32.newbyteorder(self.endianness))
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
            self._Q[grid_index] = np.zeros(grid_size + [5], dtype = self.scalarType)
        return None

    def Export(self, filename, Q_header = None):
        f = open(filename, "wb")
        record_size = self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, record_size))
        f.write(struct.pack(self.integerTypeStr, self.nGrids))
        f.write(struct.pack(self.offsetTypeStr, record_size))
        record_size = 3 * self.nGrids * self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, record_size))
        for iGrid in range(0, self.nGrids):
            gridSize = np.array(self.GetSize(iGrid), dtype = self.integerType)
            f.write(gridSize.tostring())
        f.write(struct.pack(self.offsetTypeStr, record_size))
        if Q_header is None:
            Q_header = np.zeros(4, dtype = self.scalarType)
        for iGrid in range(0, self.nGrids):
            record_size = 4 * self.scalarType.itemsize
            f.write(struct.pack(self.offsetTypeStr, record_size))
            f.write(Q_header.tostring())
            f.write(struct.pack(self.offsetTypeStr, record_size))
            record_size = 5 * self.scalarType.itemsize * np.product(self.GetSize(iGrid))
            f.write(struct.pack(self.offsetTypeStr, record_size))
            f.write(self._Q[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.offsetTypeStr, record_size))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType))
        f.seek(self.offsetType.itemsize, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(3 * self.offsetType.itemsize + 4 * self.scalarType.itemsize, 1)
            dtype = 5 * np.product(self.GetSize(iGrid)) * self.scalarType
            self._Q[iGrid][:,:,:,:] = np.reshape(np.fromstring(f.read(dtype.itemsize), dtype = dtype), self.GetSize(iGrid) + [5], order = 'F')
            f.seek(self.offsetType.itemsize, 1)
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
            f.seek(4 * self.offsetBytes + (4 + 5 * np.product(self.GetSize(iGrid))) * self.scalarBytes, 1)
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
                f.seek(5 * np.product(grid_size[iGrid]) * self.scalarBytes, 1)
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
                            self._Q[0][:, j - start_indices[1], k - start_indices[2], l] = np.array(struct.unpack(self.endianness + str((end_indices[0] - start_indices[0] + 1)) + self.scalarType, f.read(self.scalarBytes * (end_indices[0] - start_indices[0] + 1))), dtype = self.scalarType.newbyteorder(self.endianness))
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
            self._F[grid_index] = np.zeros(grid_size + [self._nComponents], dtype = self.scalarType)
        return None

    def Export(self, filename):
        f = open(filename, "wb")
        record_size = self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, record_size))
        f.write(struct.pack(self.integerTypeStr, self.nGrids))
        f.write(struct.pack(self.offsetTypeStr, record_size))
        record_size = 4 * self.nGrids * self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, record_size))
        for iGrid in range(0, self.nGrids):
            gridSize = np.array(self.GetSize(iGrid), dtype = self.integerType)
            f.write(gridSize.tostring())
            f.write(struct.pack(self.integerTypeStr, self._nComponents))
        f.write(struct.pack(self.offsetTypeStr, record_size))
        for iGrid in range(0, self.nGrids):
            record_size = self._nComponents * self.scalarType.itemsize * np.product(self.GetSize(iGrid))
            f.write(struct.pack(self.offsetTypeStr, record_size))
            f.write(self._F[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.offsetTypeStr, record_size))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(0, self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType))
            if iGrid == 0:
                self._nComponents, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
            else:
                number_of_components, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
                if number_of_components != self._nComponents:
                    raise IOError("%s: Invalid PLOT3D function file: number of components in block %i (= %i) differs from %i" % (filename, iGrid + 1, number_of_components, self._nComponents))                
            self.Update(iGrid)
        f.seek(self.offsetType.itemsize, 1)
        for iGrid in range(0, self.nGrids):
            f.seek(self.offsetType.itemsize, 1)
            dtype = self._nComponents * np.product(self.GetSize(iGrid)) * self.scalarType
            self._F[iGrid][:,:,:,:] = np.reshape(np.fromstring(f.read(dtype.itemsize), dtype = dtype), self.GetSize(iGrid) + [self._nComponents], order = 'F')
            f.seek(self.offsetType.itemsize, 1)
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
            f.seek(2 * self.offsetBytes + (self._nComponents * np.product(self.GetSize(iGrid))) * self.scalarType, 1)
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
                f.seek(self._nComponents * np.product(grid_size[iGrid]) * self.scalarBytes, 1)
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
                            self._F[0][:, j - start_indices[1], k - start_indices[2], l] = np.array(struct.unpack(self.endianness + str((end_indices[0] - start_indices[0] + 1)) + self.scalarType, f.read(self.scalarBytes * (end_indices[0] - start_indices[0] + 1))), dtype = self.scalarType.newbyteorder(self.endianness))
                            f.seek(self.scalarBytes * (grid_size[iGrid][0] - (end_indices[0] + 1)), 1)
                        f.seek(self.scalarBytes * grid_size[iGrid][0] * (grid_size[iGrid][1] - (end_indices[1] + 1)), 1)
                    f.seek(self.scalarBytes * grid_size[iGrid][0] * grid_size[iGrid][1] * (grid_size[iGrid][2] - (end_indices[2] + 1)), 1)
                if show_progress is True:
                    pbar.finish()
            f.seek(self.offsetBytes, 1)
        f.close()

        return None    
