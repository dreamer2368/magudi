#!/usr/bin/env python
import sys
import random
import numpy as np
import matplotlib.pyplot as plt

import PLOT3D

class InstabilityMode:

    def __init__(self):

        self.EPSILON = np.finfo(float).eps

        self.azimuthalWavenumber = 0
        self.temporalAngularFrequency = 0.

        self.momentumThickness = 1.
        self.machNumber = 0.
        self.nozzleLipRadius = 0.5
        self.ratioOfSpecificHeats = 1.4

    def update(self):

        self._temperatureRatio = 1. / (1. + 0.5 * (self.ratioOfSpecificHeats - 1.) * self.machNumber ** 2)
        self._streamwiseVelocityAtExit = self.machNumber * np.sqrt(self._temperatureRatio)

    def setDomain(self, radialCoordinate, rLower, rMatch, rUpper):

        self._r = radialCoordinate

        iLower = np.argmin(np.abs(self._r - rLower))
        iMatch = np.argmin(np.abs(self._r - rMatch))
        iUpper = np.argmin(np.abs(self._r - rUpper))

        self._nStepsInner = iMatch - iLower
        self._nStepsOuter = iUpper - iMatch
        
        self._rLower = self._r[iLower]
        self._rMatch = self._r[iMatch]
        self._rUpper = self._r[iUpper]

    def meanFlow(self, r):

        normalizedExitVelocity = 0.5 * (1. + np.tanh(0.25 / self.momentumThickness * (self.nozzleLipRadius / (r + self.EPSILON) - r / self.nozzleLipRadius)))
        derivativeOfNormalizedExitVelocity = -0.125 / self.momentumThickness * (1. - np.tanh(0.25 / self.momentumThickness * (self.nozzleLipRadius / (r + self.EPSILON) - r / self.nozzleLipRadius)) ** 2) * (self.nozzleLipRadius / (r + self.EPSILON) ** 2 + 1. / self.nozzleLipRadius)

        # u, \rho, a, du/dr, d\rho/dr
        meanFlow = np.array([np.empty_like(r) for i in range(5)])

        meanFlow[0] = self.machNumber * np.sqrt(self._temperatureRatio) * normalizedExitVelocity
        meanFlow[1] = 1. / self._temperatureRatio / (0.5 * (self.ratioOfSpecificHeats - 1.) * normalizedExitVelocity * (1. - normalizedExitVelocity) * self.machNumber ** 2 + normalizedExitVelocity + (1. - normalizedExitVelocity) / self._temperatureRatio)
        meanFlow[2] = np.sqrt(1. / meanFlow[1])
        meanFlow[3] = self.machNumber * np.sqrt(self._temperatureRatio) * derivativeOfNormalizedExitVelocity
        meanFlow[4] = -self._temperatureRatio * meanFlow[1] ** 2 * derivativeOfNormalizedExitVelocity * (0.5 * (self.ratioOfSpecificHeats - 1.) * (1. - 2. * normalizedExitVelocity) * self.machNumber ** 2 + 1. - 1. / self._temperatureRatio)

        return meanFlow
        
    def rhsCompressibleRayleigh(self, r, pressureEigenfunction, streamwiseWavenumber):
        u, rho, a, u_r, rho_r = self.meanFlow(r)
        return np.array([pressureEigenfunction[1], 
                         - (1. / r - rho_r / rho + 2. * streamwiseWavenumber / (self.temporalAngularFrequency - streamwiseWavenumber * u) * u_r) * pressureEigenfunction[1] - ((self.temporalAngularFrequency - streamwiseWavenumber * u) ** 2 / a ** 2 - self.azimuthalWavenumber ** 2 / r ** 2 - streamwiseWavenumber ** 2) * pressureEigenfunction[0]])

    def innerSolution(self, r, streamwiseWavenumber):
        import scipy.special
        eta = np.sqrt((self.temporalAngularFrequency - streamwiseWavenumber * self._streamwiseVelocityAtExit) ** 2 / self._temperatureRatio - streamwiseWavenumber ** 2)
        return np.array([scipy.special.jn(self.azimuthalWavenumber, eta * r), 
                         eta * 0.5 * (scipy.special.jn(self.azimuthalWavenumber - 1, eta * r) - scipy.special.jn(self.azimuthalWavenumber + 1, eta * r))])

    def outerSolution(self, r, streamwiseWavenumber):
        import scipy.special
        eta = np.sqrt(self.temporalAngularFrequency ** 2 - streamwiseWavenumber ** 2)
        return np.array([scipy.special.hankel1(self.azimuthalWavenumber, np.sign(eta.imag) * eta * r), 
                         np.sign(eta.imag) * eta * 0.5 * (scipy.special.hankel1(self.azimuthalWavenumber - 1, np.sign(eta.imag) * eta * r) - scipy.special.hankel1(self.azimuthalWavenumber + 1, np.sign(eta.imag) * eta * r))])

    def derivativeMatchError(self, streamwiseWavenumber):

        from scipy.integrate import ode

        integrator = ode(self.rhsCompressibleRayleigh)
        integrator.set_integrator('zvode', method = 'bdf', order = 15).set_f_params(streamwiseWavenumber)

        # Inner solution
        integrator.set_initial_value(self.innerSolution(self._rLower, streamwiseWavenumber), self._rLower)
        dr = (self._rMatch - self._rLower) / self._nStepsInner
        while integrator.successful() and integrator.t < self._rMatch - 0.5 * dr:
            integrator.integrate(integrator.t + dr)
        innerSolution = integrator.y

        # Outer solution
        integrator.set_initial_value(self.outerSolution(self._rUpper, streamwiseWavenumber), self._rUpper)
        dr = (self._rUpper - self._rMatch) / self._nStepsOuter
        while integrator.successful() and integrator.t > self._rMatch + 0.5 * dr:
            integrator.integrate(integrator.t - dr)
        outerSolution = integrator.y

        if abs(outerSolution[0] / innerSolution[0]) < 1e-6:
            return innerSolution[1] - outerSolution[1] * innerSolution[0] / outerSolution[0]
        return innerSolution[1] * outerSolution[0] / innerSolution[0] - outerSolution[1]

    def findEigenvalue(self, initialGuess, maxIterations = 200, tolerance = 1e-10):
        x = initialGuess
        y = self.derivativeMatchError(x)
        it = 1
        while abs(y) > tolerance and it <= maxIterations: # Newton-Raphson
            dx = random.uniform(1e-10, 1e-8) + 1j * random.uniform(1e-10, 1e-8) # step size for computing a derivative
            # print "Newton-Raphson: it = ", it, ", x = ", x, ", y = ", y, ", dx = ", dx
            x = x - dx * y / (self.derivativeMatchError(x + dx) - y)
            y = self.derivativeMatchError(x)
            it = it + 1
        if it >= maxIterations:
            print "Newton-Raphson failed: initial guess = ", initialGuess, ", current guess = ", x, ", current function = ", abs(y)
        self.streamwiseWavenumber = x.real - 1.j * abs(x.imag)

    def findEigenfunction(self, tolerance = 1e-8):

        from scipy.integrate import ode
        from scipy.interpolate import interp1d

        integrator = ode(self.rhsCompressibleRayleigh)
        integrator.set_integrator('zvode', method = 'bdf', order = 15).set_f_params(self.streamwiseWavenumber)

        pressureEigenfunction = np.array([np.empty_like(self._r) for i in range(2)], dtype = 'complex128')

        iLower = np.argmin(np.abs(self._r - self._rLower))
        pressureEigenfunction[:,:iLower] = self.innerSolution(self._r[:iLower], self.streamwiseWavenumber)

        i = iLower
        integrator.set_initial_value(self.innerSolution(self._rLower, self.streamwiseWavenumber), self._rLower)
        dr = (self._rMatch - self._rLower) / self._nStepsInner
        while integrator.successful() and integrator.t < self._rMatch - 0.5 * dr:
            pressureEigenfunction[:,i] = integrator.y
            integrator.integrate(integrator.t + dr)
            i += 1
        pressureEigenfunction[:,i] = integrator.y

        iUpper = np.argmin(np.abs(self._r - self._rUpper))
        pressureEigenfunction[:,iUpper:] = self.outerSolution(self._r[iUpper:], self.streamwiseWavenumber)

        i = iUpper
        integrator.set_initial_value(self.outerSolution(self._rUpper, self.streamwiseWavenumber), self._rUpper)
        dr = (self._rUpper - self._rMatch) / self._nStepsOuter
        while integrator.successful() and integrator.t > self._rMatch + 0.5 * dr:
            pressureEigenfunction[:,i] = integrator.y
            integrator.integrate(integrator.t - dr)
            i -= 1
        outerSolution = integrator.y

        scalingFactor = integrator.y[0] / pressureEigenfunction[0,i]
        scalingFactorInverse = 1. / scalingFactor
        if abs(scalingFactor) < 1e-6:
            pressureEigenfunction[:,i+1:] *= scalingFactorInverse
        else:
            pressureEigenfunction[:,:i+1] *= scalingFactor
            
        assert abs(pressureEigenfunction[1,i] - integrator.y[1]) < tolerance

        self.pressureEigenfunctionDerivative = interp1d(self._r, pressureEigenfunction[1,:], kind = 'linear')
        self.pressureEigenfunction = interp1d(self._r, pressureEigenfunction[0,:], kind = 'linear')

    def getEigenmode(self, radialCoordinate):
        u, rho, a, u_r, rho_r = self.meanFlow(radialCoordinate)
        p_hat = self.pressureEigenfunction(radialCoordinate)
        p_hat_r = self.pressureEigenfunctionDerivative(radialCoordinate)
        omega = self.temporalAngularFrequency - self.streamwiseWavenumber * u
        u_hat_r = 1. / (1.j * rho * omega) * p_hat_r
        u_hat_x = self.streamwiseWavenumber * p_hat / (rho * omega) + u_hat_r / (1.j * omega) * u_r
        rho_hat = p_hat / a ** 2 + u_hat_r / (1.j * omega) * rho_r
        return np.array([rho_hat,
                         rho * u_hat_r,
                         self.azimuthalWavenumber * p_hat / (radialCoordinate * omega + self.EPSILON),
                         rho_hat * u + rho * u_hat_x,
                         p_hat / (self.ratioOfSpecificHeats - 1.) + 0.5 * rho_hat * u ** 2 + rho * u * u_hat_x])

if __name__ == '__main__':

    outputPrefix = 'OSUMach1.3'
    gridFile = outputPrefix + '.xyz'

    machNumber = 1.3
    momentumThickness = 0.02

    mode = InstabilityMode()
    mode.machNumber = machNumber
    mode.momentumThickness = momentumThickness
    mode.update()
    streamwiseVelocityAtExit = mode.meanFlow(mode.EPSILON)[0]

    azimuthalWavenumbers = np.arange(1, 6)
    strouhalNumbers = np.array([0.43, 0.51, 0.61, 0.69, 0.74, 0.88])
    streamwiseWavenumbers = np.empty([azimuthalWavenumbers.size, strouhalNumbers.size], dtype = 'complex128')
    phaseAngles = np.empty([azimuthalWavenumbers.size, strouhalNumbers.size, 2])

    rLower = 0.025
    rMatch = 0.25
    rUpper = 3.
    radialCoordinate = np.linspace(0., rLower, 5001)
    radialCoordinate = np.append(radialCoordinate[:-1], np.linspace(rLower, rMatch, 5001))
    radialCoordinate = np.append(radialCoordinate[:-1], np.linspace(rMatch, rUpper, 5001))
    radialCoordinate = np.append(radialCoordinate[:-1], np.linspace(rUpper, 15., 5001))

    showProgress = True
    try:
        import progressbar
    except:
        showProgress = False

    if showProgress:
        progressBar = progressbar.ProgressBar(maxval = streamwiseWavenumbers.size)
        progressBar.start()

    modes = [[InstabilityMode() for azimuthalWavenumber in azimuthalWavenumbers] for strouhalNumber in strouhalNumbers]

    initialGuess = 1. - 1.j
    for j, strouhalNumber in enumerate(strouhalNumbers):
        guess = initialGuess
        for i, azimuthalWavenumber in enumerate(azimuthalWavenumbers):
            modes[j][i].machNumber = machNumber
            modes[j][i].momentumThickness = momentumThickness
            modes[j][i].update()
            modes[j][i].setDomain(radialCoordinate, rLower, rMatch, rUpper)
            modes[j][i].azimuthalWavenumber = azimuthalWavenumber
            modes[j][i].temporalAngularFrequency = 2. * np.pi * streamwiseVelocityAtExit * strouhalNumber
            modes[j][i].findEigenvalue(guess)
            modes[j][i].findEigenfunction()
            phaseAngles[i,j,0] = 2. * np.pi * random.random()
            phaseAngles[i,j,1] = 2. * np.pi * random.random()
            guess = streamwiseWavenumbers[i,j] = modes[j][i].streamwiseWavenumber
            if i == 0:
                initialGuess = modes[j][i].streamwiseWavenumber
            if showProgress:
                progressBar.update(i + azimuthalWavenumbers.size * j)

    if showProgress:
        progressBar.finish()

    grid = PLOT3D.Grid()
    grid.hasIBLANK = False
    grid.ImportSkeleton(gridFile)
    nGrids = grid.nGrids
    assert nGrids == 5

    for k, strouhalNumber in enumerate(strouhalNumbers):

        excitationPatchGrid = PLOT3D.Grid()
        excitationPatchGrid.nGrids = nGrids

        solnReal = PLOT3D.Solution()
        solnReal.nGrids = nGrids
        solnImag = PLOT3D.Solution()
        solnImag.nGrids = nGrids

        for i in range(nGrids):

            grid.ImportSubarray(gridFile, i, [0, 0, 0], [0, 0, 32])
            axialCoordinate = grid.X[0][0,0,:,2]
            grid.ImportSubarray(gridFile, i, [0, 0, 0], [-1, -1, 0])
            radialCoordinate = np.sqrt(grid.X[0][:,:,0,0] ** 2 + grid.X[0][:,:,0,1] ** 2)
            polarAngle = np.arctan2(grid.X[0][:,:,0,1], grid.X[0][:,:,0,0])
            excitationPatchGrid.SetSize(i, [ radialCoordinate.shape[0], radialCoordinate.shape[1], axialCoordinate.size ])
            solnReal.SetSize(i, excitationPatchGrid.GetSize(i))
            solnImag.SetSize(i, excitationPatchGrid.GetSize(i))
            solnReal.auxiliaryData[i][1] = solnImag.auxiliaryData[i][1] = modes[k][0].temporalAngularFrequency
            for k_ in range(axialCoordinate.size):
                excitationPatchGrid.X[i][:,:,k_,0:2] = grid.X[0][:,:,0,0:2]
                excitationPatchGrid.X[i][:,:,k_,2] = axialCoordinate[k_]
        
            for j, azimuthalWavenumber in enumerate(azimuthalWavenumbers):
                Q = modes[k][j].getEigenmode(radialCoordinate)
                for l in range(Q.shape[0]):
                    if l == 2:
                        Q[l] = Q[l] * (np.exp(1.j * azimuthalWavenumber * polarAngle + phaseAngles[j,k,0]) - np.exp(-1.j * azimuthalWavenumber * polarAngle + phaseAngles[j,k,1]))
                    else:
                        Q[l] = Q[l] * (np.exp(1.j * azimuthalWavenumber * polarAngle + phaseAngles[j,k,0]) + np.exp(-1.j * azimuthalWavenumber * polarAngle + phaseAngles[j,k,1]))
                U = Q[1] * np.cos(polarAngle) - Q[2] * np.sin(polarAngle)
                V = Q[1] * np.sin(polarAngle) + Q[2] * np.cos(polarAngle)
                Q[1] = U
                Q[2] = V
                for l in range(Q.shape[0]):
                    for k_ in range(axialCoordinate.size):
                        solnReal.Q[i][:,:,k_,l] += np.real(Q[l] * np.exp(1.j * modes[k][j].streamwiseWavenumber * axialCoordinate[k_]))
                        solnImag.Q[i][:,:,k_,l] += np.imag(Q[l] * np.exp(1.j * modes[k][j].streamwiseWavenumber * axialCoordinate[k_]))

        if k == 0:
            excitationPatchGrid.Export(outputPrefix + '.eigenmodes.xyz')
        solnReal.Export(outputPrefix + '-%02d.eigenmode_real.q' % (k + 1))
        solnImag.Export(outputPrefix + '-%02d.eigenmode_imag.q' % (k + 1))
