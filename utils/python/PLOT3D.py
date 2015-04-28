#!/usr/bin/env python
import numpy as np
import struct

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
        allowedValues = ['=', '<', '>']
        if value not in allowedValues:
            raise ValueError("'%s' is not a valid endianness indicator (allowed values are: %s)" % (value, ", ".join(map(str, allowedValues))))
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
        allowedValues = np.sctypes['float'] + np.sctypes['complex']
        if value not in allowedValues:
            raise ValueError("'%s' is not a valid format string for a scalar type (allowed values are: %s)" % (value, ", ".join(map(str, allowedValues))))
        self._scalarType = np.dtype(value).newbyteorder(self._endianness)

    @property
    def integerType(self):
        return self._integerType

    @integerType.setter
    def integerType(self, value):
        allowedValues = np.sctypes['int'] + np.sctypes['uint']
        if value not in allowedValues:
            raise ValueError("'%s' is not a valid format string for a integer type (allowed values are: %s)" % (value, ", ".join(map(str, allowedValues))))
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
        allowedValues = np.sctypes['uint']
        if value not in allowedValues:
            raise ValueError("'%s' is not a valid format string for an offset type (allowed values are: %s)" % (value, ", ".join(map(str, allowedValues))))
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
        self._size = [ [0, 0, 0] for iGrid in range(self._nGrids) ]
        self.Update()

    @property
    def size(self):
        return self._size

    def GetSize(self, gridIndex = None):
        if gridIndex is None:
            return self._size
        if not isinstance(gridIndex, (int, long)):
            raise TypeError("Grid index must be an integer")
        if gridIndex < 0 or gridIndex >= self._nGrids:
            raise IndexError("Grid index out of range (expected a number in [%i,%i])" % (0, self._nGrids - 1))
        return self._size[gridIndex]

    def SetSize(self, gridIndex, size, allocate = True):
        if not isinstance(gridIndex, (int, long)):
            raise TypeError("Grid index must be an integer")
        if gridIndex < 0 or gridIndex >= self._nGrids:
            raise IndexError("Grid index out of range (expected a number in [%i,%i])" % (0, self._nGrids - 1))
        assert len(size) == 3
        for i in range(3):
            self._size[gridIndex][i] = size[i]
        if allocate is True:
            self.Update(gridIndex)
        return None

    def GetNumberOfPoints(self):
        nPoints = 0
        for iGrid in range(self._nGrids):
            nPoints = nPoints + np.product(self._size[iGrid])
        return nPoints

    def CopyFrom(self, obj):
        self._endianness = obj.endianness
        self._integerType = obj.integerType
        self._scalarType = obj.scalarType
        self._offsetType = obj.offsetType
        self._integerTypeStr = obj.integerTypeStr
        self._offsetTypeStr = obj.offsetTypeStr
        if obj.nGrids > 0:
            self.nGrids = obj.nGrids
            for iGrid in range(self._nGrids):
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

    def Update(self, gridIndex = None):
        if gridIndex is None:
            self._X = [ None for iGrid in range(self.nGrids) ]
            self._IBLANK = [ None for iGrid in range(self.nGrids) ]
        else:
            gridSize = self.GetSize(gridIndex)
            self._X[gridIndex] = np.zeros(gridSize + [3], dtype = self.scalarType)
            self._IBLANK[gridIndex] = np.empty(gridSize, dtype = self.integerType)
            self._IBLANK[gridIndex].fill(1)
        return None

    def Export(self, filename):
        f = open(filename, "wb")
        recordSize = self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        f.write(struct.pack(self.integerTypeStr, self.nGrids))
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        recordSize = 3 * self.nGrids * self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        for iGrid in range(self.nGrids):
            gridSize = np.array(self.GetSize(iGrid), dtype = self.integerType)
            f.write(gridSize.tostring())
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        for iGrid in range(self.nGrids):
            recordSize = (self.integerType.itemsize + 3 * self.scalarType.itemsize) * np.product(self.GetSize(iGrid))
            f.write(struct.pack(self.offsetTypeStr, recordSize))
            f.write(self._X[iGrid].tostring(order = 'F'))
            f.write(self._IBLANK[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.offsetTypeStr, recordSize))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType))
        f.seek(self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
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
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType), allocate = False)
        f.close()
        return None

    def ImportSubarray(self, filename, gridIndex, startIndices, endIndices, showProgress = False):

        try:
            if showProgress:
                from progressbar import ProgressBar
        except:
            showProgress = False

        assert len(startIndices) == 3 and len(endIndices) == 3
        self.ImportSkeleton(filename)
        assert gridIndex >= 0 and gridIndex < self.nGrids

        for i in range(3):
            if startIndices[i] < 0:
                startIndices[i] += self._size[gridIndex][i]
            if endIndices[i] < 0:
                endIndices[i] += self._size[gridIndex][i]
            assert startIndices[i] >= 0 and endIndices[i] < self._size[gridIndex][i] and endIndices[i] >= startIndices[i]

        gridSize = np.copy(self.GetSize())
        nGrids = len(gridSize)

        self.nGrids = 1
        self.SetSize(0, [ endIndices[i] - startIndices[i] + 1 for i in range(3) ])

        f = open(filename, "rb")
        f.seek(4 * self.offsetType.itemsize + (3 * nGrids + 1) * self.integerType.itemsize)
        for iGrid in range(nGrids):
            f.seek(self.offsetType.itemsize, 1)
            if iGrid != gridIndex:
                f.seek(3 * np.product(gridSize[iGrid]) * self.scalarType.itemsize, 1)
                if self._hasIBLANK is True:
                    f.seek(np.product(gridSize[iGrid]) * self.integerType.itemsize, 1)
            else:
                
                if showProgress is True:
                    if self._hasIBLANK is True:
                        progressBar = ProgressBar(maxval = 4 * (endIndices[2] - startIndices[2] + 1) * (endIndices[1] - startIndices[1] + 1))
                    else:
                        progressBar = ProgressBar(maxval = 3 * (endIndices[2] - startIndices[2] + 1) * (endIndices[1] - startIndices[1] + 1))
                    progressBar.start()

                for i in range(3):
                    f.seek(gridSize[iGrid][0] * gridSize[iGrid][1] * startIndices[2] * self.scalarType.itemsize, 1)
                    for k in range(startIndices[2], endIndices[2] + 1):
                        f.seek(gridSize[iGrid][0] * startIndices[1] * self.scalarType.itemsize, 1)
                        for j in range(startIndices[1], endIndices[1] + 1):
                            if showProgress is True:
                                progressBar.update(j - startIndices[1] + 1 + (endIndices[1] - startIndices[1] + 1) * (k - startIndices[2] + (endIndices[2] - startIndices[2] + 1) * i))
                            f.seek(startIndices[0] * self.scalarType.itemsize, 1)
                            self._X[0][:, j - startIndices[1], k - startIndices[2], i] = np.fromstring(f.read((endIndices[0] - startIndices[0] + 1) * self.scalarType.itemsize), dtype = self.scalarType)
                            f.seek((gridSize[iGrid][0] - (endIndices[0] + 1)) * self.scalarType.itemsize, 1)
                        f.seek(gridSize[iGrid][0] * (gridSize[iGrid][1] - (endIndices[1] + 1)) * self.scalarType.itemsize, 1)
                    f.seek(gridSize[iGrid][0] * gridSize[iGrid][1] * (gridSize[iGrid][2] - (endIndices[2] + 1)) * self.scalarType.itemsize, 1)

                if self._hasIBLANK is True:
                    f.seek(gridSize[iGrid][0] * gridSize[iGrid][1] * startIndices[2] * self.integerType.itemsize, 1)
                    for k in range(startIndices[2], endIndices[2] + 1):
                        f.seek(gridSize[iGrid][0] * startIndices[1] * self.integerType.itemsize, 1)
                        for j in range(startIndices[1], endIndices[1] + 1):
                            if showProgress is True:
                                progressBar.update(j - startIndices[1] + 1 + (endIndices[1] - startIndices[1] + 1) * (k - startIndices[2] + (endIndices[2] - startIndices[2] + 1) * 3))
                            f.seek(startIndices[0] * self.integerType.itemsize, 1)
                            self._IBLANK[0][:, j - startIndices[1], k - startIndices[2]] = np.fromstring(f.read((endIndices[0] - startIndices[0] + 1) * self.integerType.itemsize), dtype = self.integerType)
                            f.seek((gridSize[iGrid][0] - (endIndices[0] + 1)) * self.integerType.itemsize, 1)
                        f.seek(gridSize[iGrid][0] * (gridSize[iGrid][1] - (endIndices[1] + 1)) * self.integerType.itemsize, 1)
                    f.seek(gridSize[iGrid][0] * gridSize[iGrid][1] * (gridSize[iGrid][2] - (endIndices[2] + 1)) * self.integerType.itemsize, 1)

                if showProgress is True:
                    progressBar.finish()

            f.seek(self.offsetType.itemsize, 1)
        f.close()

        return None

    def CopyFrom(self, obj):
        super(Grid, self).CopyFrom(obj)
        try:
            self._hasIBLANK = obj._hasIBLANK
        except:
            pass

class Solution(MultiBlockObject):

    def __init__(self):
        super(Solution, self).__init__()
        self._Q = []
        self._auxiliaryData = []

    @property
    def Q(self):
        return self._Q

    @property
    def auxiliaryData(self):
        return self._auxiliaryData

    def Update(self, gridIndex = None):
        if gridIndex is None:
            self._Q = [ None ] * self.nGrids
            self._auxiliaryData = [ np.zeros(4, dtype = self.scalarType) ] * self.nGrids
        else:
            gridSize = self.GetSize(gridIndex)
            self._Q[gridIndex] = np.zeros(gridSize + [5], dtype = self.scalarType)
        return None

    def Export(self, filename):
        f = open(filename, "wb")
        recordSize = self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        f.write(struct.pack(self.integerTypeStr, self.nGrids))
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        recordSize = 3 * self.nGrids * self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        for iGrid in range(self.nGrids):
            gridSize = np.array(self.GetSize(iGrid), dtype = self.integerType)
            f.write(gridSize.tostring())
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        for iGrid in range(self.nGrids):
            recordSize = 4 * self.scalarType.itemsize
            f.write(struct.pack(self.offsetTypeStr, recordSize))
            f.write(self._auxiliaryData[iGrid].tostring())
            f.write(struct.pack(self.offsetTypeStr, recordSize))
            recordSize = 5 * self.scalarType.itemsize * np.product(self.GetSize(iGrid))
            f.write(struct.pack(self.offsetTypeStr, recordSize))
            f.write(self._Q[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.offsetTypeStr, recordSize))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType))
        f.seek(self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            f.seek(self.offsetType.itemsize, 1)
            self._auxiliaryData[iGrid] = np.fromstring(f.read(4 * self.scalarType.itemsize), dtype = self.scalarType)
            f.seek(2 * self.offsetType.itemsize, 1)
            dtype = 5 * np.product(self.GetSize(iGrid)) * self.scalarType
            self._Q[iGrid][:,:,:,:] = np.reshape(np.fromstring(f.read(dtype.itemsize), dtype = dtype), self.GetSize(iGrid) + [5], order = 'F')
            f.seek(self.offsetType.itemsize, 1)
        f.close()
        return None

    def ImportSkeleton(self, filename):
        f = open(filename, "rb")
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType), allocate = False)
        f.seek(self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            f.seek(self.offsetType.itemsize, 1)
            self._auxiliaryData[iGrid] = np.fromstring(f.read(4 * self.scalarType.itemsize), dtype = self.scalarType)
            dtype = 5 * np.product(self.GetSize(iGrid)) * self.scalarType
            f.seek(3 * self.offsetType.itemsize + dtype.itemsize, 1)
        f.close()
        return None

    def ImportSubarray(self, filename, gridIndex, startIndices, endIndices, showProgress = False):

        try:
            if showProgress:
                from progressbar import ProgressBar
        except:
            showProgress = False

        assert len(startIndices) == 3 and len(endIndices) == 3
        self.ImportSkeleton(filename)
        assert gridIndex >= 0 and gridIndex < self.nGrids

        for i in range(3):
            if startIndices[i] < 0:
                startIndices[i] += self._size[gridIndex][i]
            if endIndices[i] < 0:
                endIndices[i] += self._size[gridIndex][i]
            assert startIndices[i] >= 0 and endIndices[i] < self._size[gridIndex][i] and endIndices[i] >= startIndices[i]

        gridSize = np.copy(self.GetSize())
        nGrids = gridSize
 
        self.nGrids = 1
        self.SetSize(0, [ endIndices[i] - startIndices[i] + 1 for i in range(3) ])

        f = open(filename, "rb")
        f.seek(4 * self.offsetType.itemsize + (3 * nGrids + 1) * self.integerType.itemsize)
        for iGrid in range(nGrids):
            if iGrid != gridIndex:
                f.seek(3 * self.offsetType.itemsize + 4 * self.scalarType.itemsize + 5 * np.product(gridSize[iGrid]) * self.scalarType.itemsize, 1)
            else:

                f.seek(self.offsetType.itemsize, 1)
                self._auxiliaryData[0] = np.fromstring(f.read(4 * self.scalarType.itemsize), dtype = self.scalarType)
                f.seek(2 * self.offsetType.itemsize, 1)
                
                if showProgress is True:
                    progressBar = ProgressBar(maxval = 5 * (endIndices[2] - startIndices[2] + 1) * (endIndices[1] - startIndices[1] + 1))
                    progressBar.start()

                for i in range(5):
                    f.seek(gridSize[iGrid][0] * gridSize[iGrid][1] * startIndices[2] * self.scalarType.itemsize, 1)
                    for k in range(startIndices[2], endIndices[2] + 1):
                        f.seek(gridSize[iGrid][0] * startIndices[1] * self.scalarType.itemsize, 1)
                        for j in range(startIndices[1], endIndices[1] + 1):
                            if showProgress is True:
                                progressBar.update(j - startIndices[1] + 1 + (endIndices[1] - startIndices[1] + 1) * (k - startIndices[2] + (endIndices[2] - startIndices[2] + 1) * i))
                            f.seek(startIndices[0] * self.scalarType.itemsize, 1)
                            self._Q[0][:, j - startIndices[1], k - startIndices[2], i] = np.fromstring(f.read((endIndices[0] - startIndices[0] + 1) * self.scalarType.itemsize), dtype = self.scalarType)
                            f.seek((gridSize[iGrid][0] - (endIndices[0] + 1)) * self.scalarType.itemsize, 1)
                        f.seek(gridSize[iGrid][0] * (gridSize[iGrid][1] - (endIndices[1] + 1)) * self.scalarType.itemsize, 1)
                    f.seek(gridSize[iGrid][0] * gridSize[iGrid][1] * (gridSize[iGrid][2] - (endIndices[2] + 1)) * self.scalarType.itemsize, 1)

                if showProgress is True:
                    progressBar.finish()

            f.seek(self.offsetType.itemsize, 1)
        f.close()

        return None

    def ImportAverage(self, filename, axis, showProgress = False):

        try:
            if showProgress:
                from progressbar import ProgressBar
        except:
            showProgress = False

        assert axis in [0, 1, 2]
        self.ImportSkeleton(filename)

        gridSize = np.copy(self.GetSize())
        for i in range(self.nGrids):
            gridSize_ = self.GetSize(i)
            gridSize_[axis] = 1
            self.SetSize(i, gridSize_)
 
        f = open(filename, "rb")
        f.seek(4 * self.offsetType.itemsize + (3 * self.nGrids + 1) * self.integerType.itemsize)

        if showProgress is True:
            progressBar = ProgressBar(maxval = 5 * self.nGrids * gridSize[iGrid][1] * gridSize[iGrid][2])
            progressBar.start()

        for iGrid in range(self.nGrids):
            f.seek(self.offsetType.itemsize, 1)
            self._auxiliaryData[iGrid] = np.fromstring(f.read(4 * self.scalarType.itemsize), dtype = self.scalarType)
            f.seek(2 * self.offsetType.itemsize, 1)
            for l in range(5):
                if axis == 0:
                    for k in range(gridSize[iGrid][2]):
                        for j in range(gridSize[iGrid][1]):
                            self._Q[iGrid][0,j,k,l] = np.mean(np.fromstring(f.read(gridSize[iGrid][0] * self.scalarType.itemsize), dtype = self.scalarType))
                            if showProgress is True:
                                progressBar.update(j + gridSize[iGrid][1] * (k + gridSize[iGrid][2] * (l + 5 * iGrid)))
                elif axis == 1:
                    for k in range(gridSize[iGrid][2]):
                        for j in range(gridSize[iGrid][1]):
                            self._Q[iGrid][:,0,k,l] += np.fromstring(f.read(gridSize[iGrid][0] * self.scalarType.itemsize), dtype = self.scalarType)                    
                            if showProgress is True:
                                progressBar.update(j + gridSize[iGrid][1] * (k + gridSize[iGrid][2] * (l + 5 * iGrid)))
                    self._Q[iGrid][:,:,:,l] /= gridSize[iGrid][1]
                else:
                    for k in range(gridSize[iGrid][2]):
                        self._Q[iGrid][:,:,0,l] += np.reshape(np.fromstring(f.read(gridSize[iGrid][0] * gridSize[iGrid][1] * self.scalarType.itemsize), dtype = self.scalarType), [gridSize[iGrid][0], gridSize[iGrid][1]], order = 'F')
                        if showProgress is True:
                            progressBar.update(gridSize[iGrid][1] * (k + 1 + gridSize[iGrid][2] * (l + 5 * iGrid)) - 1)
                    self._Q[iGrid][:,:,:,l] /= gridSize[iGrid][2]
            f.seek(self.offsetType.itemsize, 1)

        if showProgress is True:
            progressBar.finish()
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

    def Update(self, gridIndex = None):
        if gridIndex is None:
            self._F = [ None for iGrid in range(self.nGrids) ]
        else:
            gridSize = self.GetSize(gridIndex)
            self._F[gridIndex] = np.zeros(gridSize + [self._nComponents], dtype = self.scalarType)
        return None

    def Export(self, filename):
        f = open(filename, "wb")
        recordSize = self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        f.write(struct.pack(self.integerTypeStr, self.nGrids))
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        recordSize = 4 * self.nGrids * self.integerType.itemsize
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        for iGrid in range(self.nGrids):
            gridSize = np.array(self.GetSize(iGrid), dtype = self.integerType)
            f.write(gridSize.tostring())
            f.write(struct.pack(self.integerTypeStr, self._nComponents))
        f.write(struct.pack(self.offsetTypeStr, recordSize))
        for iGrid in range(self.nGrids):
            recordSize = self._nComponents * self.scalarType.itemsize * np.product(self.GetSize(iGrid))
            f.write(struct.pack(self.offsetTypeStr, recordSize))
            f.write(self._F[iGrid].tostring(order = 'F'))
            f.write(struct.pack(self.offsetTypeStr, recordSize))
        f.close()
        return None

    def Import(self, filename):
        f = open(filename, "rb")
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType))
            if iGrid == 0:
                self._nComponents, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
            else:
                number_of_components, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
                if number_of_components != self._nComponents:
                    raise IOError("%s: Invalid PLOT3D function file: number of components in block %i (= %i) differs from %i" % (filename, iGrid + 1, number_of_components, self._nComponents))                
            self.Update(iGrid)
        f.seek(self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            f.seek(self.offsetType.itemsize, 1)
            dtype = self._nComponents * np.product(self.GetSize(iGrid)) * self.scalarType
            self._F[iGrid][:,:,:,:] = np.reshape(np.fromstring(f.read(dtype.itemsize), dtype = dtype), self.GetSize(iGrid) + [self._nComponents], order = 'F')
            f.seek(self.offsetType.itemsize, 1)
        f.close()
        return None

    def ImportSkeleton(self, filename):
        f = open(filename, "rb")
        f.seek(self.offsetType.itemsize)
        self.nGrids, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
        f.seek(2 * self.offsetType.itemsize, 1)
        for iGrid in range(self.nGrids):
            self.SetSize(iGrid, np.fromstring(f.read(3 * self.integerType.itemsize), dtype = self.integerType), allocate = False)
            if iGrid == 0:
                self._nComponents, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
            else:
                number_of_components, = struct.unpack(self.integerTypeStr, f.read(self.integerType.itemsize))
                if number_of_components != self._nComponents:
                    raise IOError("%s: Invalid PLOT3D function file: number of components in block %i (= %i) differs from %i" % (filename, iGrid + 1, number_of_components, self._nComponents))                
        f.close()
        return None

    def ImportSubarray(self, filename, gridIndex, startIndices, endIndices, showProgress = False):

        try:
            if showProgress:
                from progressbar import ProgressBar
        except:
            showProgress = False

        assert len(startIndices) == 3 and len(endIndices) == 3
        self.ImportSkeleton(filename)
        assert gridIndex >= 0 and gridIndex < self.nGrids

        for i in range(3):
            if startIndices[i] < 0:
                startIndices[i] += self._size[gridIndex][i]
            if endIndices[i] < 0:
                endIndices[i] += self._size[gridIndex][i]
            assert startIndices[i] >= 0 and endIndices[i] < self._size[gridIndex][i] and endIndices[i] >= startIndices[i]

        gridSize = np.copy(self.GetSize())
        nGrids = len(gridSize)
 
        self.nGrids = 1
        self.SetSize(0, [ endIndices[i] - startIndices[i] + 1 for i in range(3) ])

        f = open(filename, "rb")
        f.seek(4 * self.offsetType.itemsize + (4 * nGrids + 1) * self.integerType.itemsize)
        for iGrid in range(nGrids):
            f.seek(self.offsetType.itemsize, 1)
            if iGrid != gridIndex:
                f.seek(self._nComponents * np.product(gridSize[iGrid]) * self.scalarType.itemsize, 1)
            else:
                
                if showProgress is True:
                    progressBar = ProgressBar(maxval = self._nComponents * (endIndices[2] - startIndices[2] + 1) * (endIndices[1] - startIndices[1] + 1))
                    progressBar.start()

                for i in range(self._nComponents):
                    f.seek(gridSize[iGrid][0] * gridSize[iGrid][1] * startIndices[2] * self.scalarType.itemsize, 1)
                    for k in range(startIndices[2], endIndices[2] + 1):
                        f.seek(gridSize[iGrid][0] * startIndices[1] * self.scalarType.itemsize, 1)
                        for j in range(startIndices[1], endIndices[1] + 1):
                            if showProgress is True:
                                progressBar.update(j - startIndices[1] + 1 + (endIndices[1] - startIndices[1] + 1) * (k - startIndices[2] + (endIndices[2] - startIndices[2] + 1) * i))
                            f.seek(startIndices[0] * self.scalarType.itemsize, 1)
                            self._F[0][:, j - startIndices[1], k - startIndices[2], i] = np.fromstring(f.read((endIndices[0] - startIndices[0] + 1) * self.scalarType.itemsize), dtype = self.scalarType)
                            f.seek((gridSize[iGrid][0] - (endIndices[0] + 1)) * self.scalarType.itemsize, 1)
                        f.seek(gridSize[iGrid][0] * (gridSize[iGrid][1] - (endIndices[1] + 1)) * self.scalarType.itemsize, 1)
                    f.seek(gridSize[iGrid][0] * gridSize[iGrid][1] * (gridSize[iGrid][2] - (endIndices[2] + 1)) * self.scalarType.itemsize, 1)

                if showProgress is True:
                    progressBar.finish()

            f.seek(self.offsetType.itemsize, 1)
        f.close()

        return None
