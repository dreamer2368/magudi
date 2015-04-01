#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve

def airfoilProfile(x, chordLength, thicknessRatio):
    return 5. * thicknessRatio * chordLength * (0.2969 * np.sqrt(x / chordLength) - 0.1260 * (x / chordLength) - 0.3516 * (x / chordLength) ** 2 + 0.2843 * (x / chordLength) ** 3 - 0.1015 * (x / chordLength) ** 4)

def airfoilProfileSlope(x, chordLength, thicknessRatio):
    return 5. * thicknessRatio * (0.14845 / np.sqrt(x / chordLength) - 0.1260 - 0.7032 * (x / chordLength) + 0.8529 * (x / chordLength) ** 2 - 0.406 * (x / chordLength) ** 3)

if __name__ == '__main__':

    outputPrefix = 'NACA0012'

    gridFile = outputPrefix + '_PLOT3D.x'
    fluentCaseFile = outputPrefix + '.cas'

    chordLength = 2000.
    thicknessRatio = 0.12

    angleOfAttack = 2.

    s = np.linspace(0., 1., 400)
    xAirfoil = np.cumsum(np.tanh(3.5 * (1. - s)) + np.tanh(3.5 * s) - 1.)
    xAirfoil = (xAirfoil - xAirfoil.min()) / (xAirfoil.max() - xAirfoil.min()) * chordLength
    yAirfoil = airfoilProfile(xAirfoil, chordLength, thicknessRatio)
    xTrailingEdge = xAirfoil[-1] - 0.5 * yAirfoil[-1] / airfoilProfileSlope(xAirfoil[-1], chordLength, thicknessRatio)
    xStart = 0.5 * chordLength

    gridSize = [360, 360, 1]
    minimumPercentSpacingAtTrailingEdge = 0.02
    minimumPercentSpacingAtLeadingEdge = 0.04
    minimumWallNormalPercentSpacing = 120.
    stopHeight = 150. * chordLength
    geometricStretchingRatio = 1.02

    minimumSpacingAtTrailingEdge = chordLength * minimumPercentSpacingAtTrailingEdge / 100.
    minimumSpacingAtLeadingEdge = chordLength * minimumPercentSpacingAtLeadingEdge / 100.
    minimumWallNormalSpacing = minimumWallNormalPercentSpacing / 100. * minimumSpacingAtTrailingEdge

    glyphFile = 'gridgenScript.glf'
    f = open(glyphFile, 'w')

    print >>f, """
gg::memClear
gg::aswSet "FLUENT" -dim 2
gg::defReset
gg::tolReset
gg::dispViewReset
set cwd [file dirname [info script]]
set cfd_solver "FLUENT"
"""

    print >>f, 'gg::dbCurveBegin -type CUBIC'
    for x, y in zip(xAirfoil[1:], yAirfoil[1:]):
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (x, y, 0.)
    print >>f, 'set db_NACA1 [gg::dbCurveEnd]'

    print >>f, 'gg::dbCurveBegin -type CUBIC'
    for x, y in zip(reversed(xAirfoil[1:]), reversed(yAirfoil[1:])):
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (x, -y, 0.)
    print >>f, 'set db_NACA2 [gg::dbCurveEnd]'

    print >>f, 'gg::dbCurveBegin -type CONIC -rho 0.5'
    print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (xAirfoil[1],  yAirfoil[1], 0.)
    print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (xAirfoil[1], -yAirfoil[1], 0.)
    print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (xAirfoil[0],  yAirfoil[0], 0.)
    print >>f, 'set db_NACA3 [gg::dbCurveEnd]'

    print >>f, 'gg::dbCurveBegin -type CONIC -rho 0.5'
    print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (xAirfoil[-1],  yAirfoil[-1], 0.)
    print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (xAirfoil[-1], -yAirfoil[-1], 0.)
    print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (xTrailingEdge,  0., 0.)
    print >>f, 'set db_NACA4 [gg::dbCurveEnd]'

    print >>f, """
set con_NACA1 [gg::conOnDBEnt $db_NACA1]
set con_NACA2 [gg::conOnDBEnt $db_NACA2]
set con_NACA3 [gg::conOnDBEnt $db_NACA3]
set con_NACA4 [gg::conOnDBEnt $db_NACA4]
gg::dbDelete ALL
"""

    print >>f, 'set con_NACA5 [gg::conSplit $con_NACA1 [gg::conGetPt $con_NACA1 -x %f]]' % (xStart)
    print >>f, """
set con_NACA6 [gg::conJoin $con_NACA5 $con_NACA4]
set con_NACA7 [gg::conJoin $con_NACA6 $con_NACA2]
set con_NACA8 [gg::conJoin $con_NACA7 $con_NACA3]
set con_NACA [gg::conJoin $con_NACA1 $con_NACA8]
"""

    print >>f, 'gg::conDim $con_NACA %i' % (gridSize[0])
    
    print >>f, 'set LE_pt [gg::conGetPt $con_NACA -x 0]'
    print >>f, 'gg::conSetBreakPt $con_NACA $LE_pt'
    print >>f, 'set TE_pt [gg::conGetPt $con_NACA -y 0]'
    print >>f, 'gg::conSetBreakPt $con_NACA $TE_pt'
    print >>f, 'gg::conBreakPtSpacing $con_NACA -breakpt 1 %f' % (minimumSpacingAtTrailingEdge)
    print >>f, 'gg::conBreakPtSpacing $con_NACA -breakpt 2 %f' % (minimumSpacingAtLeadingEdge)

    print >>f, 'gg::domExtrusionBegin [list $con_NACA] -edge'
    print >>f, 'gg::domExtrusionMode HYPERBOLIC'
    print >>f, 'gg::domExtrusionAtt -s_init %f -march_plane {0 0 1} -stop_height %f -growth_geometric %f -normal_count 20 -normal_relax 1 -vol_smoothing 0.25' % (minimumWallNormalSpacing, stopHeight, geometricStretchingRatio)
    print >>f, 'gg::domExtrusionStep -result ExtResult %i' % (gridSize[1] - 1)
    print >>f, 'set dom_NACA [gg::domExtrusionEnd]'
    
    print >>f, """
gg::blkBegin -type STRUCTURED
gg::faceBegin
gg::faceAddDom $dom_NACA
gg::faceEnd
set blk_NACA [gg::blkEnd]
gg::blkSpecifyIJK $blk_NACA 1 2
"""

    print >>f, 'gg::blkTransformBegin $blk_NACA -maintain_linkage'
    print >>f, 'gg::xformRotate [list 0 0 0] [list 0 0 1] %f' % (-angleOfAttack)
    print >>f, 'gg::blkTransformEnd'
    
    print >>f, 'gg::blkTransformBegin $blk_NACA -maintain_linkage'
    print >>f, 'gg::xformScale [list 0 0 0] [list %f %f %f]' % (1. / xTrailingEdge, 1. / xTrailingEdge, 1. / xTrailingEdge)
    print >>f, 'gg::blkTransformEnd'
    
    print >>f, """
set airfoil [lindex [gg::conGetAll] 0]
gg::aswSetBC [list $airfoil] "Wall"
set freestream [lindex [gg::conGetAll] 2]
gg::aswSetBC [list $freestream] "Velocity Inlet"
set overlap [lindex [gg::conGetAll] 1]
gg::aswSetBC [list $overlap] "Unspecified"
"""
    print >>f, 'gg::domExport $dom_NACA [file join $cwd \"%s\"] -style PLOT3D -form UNFORMATTED -precision DOUBLE -endian BIG' % (gridFile)
    
    print >>f, """
gg::conDelete $overlap -force
gg::conDim $airfoil 500
gg::conDim $freestream 2000
gg::domBegin -type UNSTRUCTURED
gg::edgeBegin
gg::edgeAddCon $freestream
gg::edgeEnd
gg::edgeBegin
gg::edgeAddCon $airfoil
gg::edgeReorient
gg::edgeEnd
set dom_NACA [gg::domEnd]
gg::blkBegin -type UNSTRUCTURED
gg::faceBegin
gg::faceAddDom $dom_NACA
gg::faceEnd
set blk_NACA [gg::blkEnd]
gg::blkTransformBegin $blk_NACA -maintain_linkage
gg::xformMirror Z
gg::blkTransformEnd
gg::aswSetBC [list $airfoil] "Wall"
gg::aswSetBC [list $freestream] "Velocity Inlet"
"""
    
    print >>f, 'gg::aswExport [file join $cwd \"%s\"]' % (fluentCaseFile)

    f.close()
