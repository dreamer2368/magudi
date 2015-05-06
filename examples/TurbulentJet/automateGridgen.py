#!/usr/bin/env python
import numpy as np
import math
import subprocess

def drawCircularArc(f, label, p1, p2, c):
    print >>f, 'gg::conBegin'
    print >>f, 'gg::segBegin -type CIRCULAR_ARC'
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p1[0], p1[1], p1[2])
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p2[0], p2[1], p2[2])
    print >>f, 'gg::segAddControlPt -alternate CENTER {%f %f %f}' % (c[0], c[1], c[2])
    print >>f, 'gg::segEnd'
    print >>f, 'set %s [gg::conEnd]' % (label)
    return None

def drawLine3D(f, label, p1, p2):
    print >>f, 'gg::conBegin'
    print >>f, 'gg::segBegin -type 3D_LINE'
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p1[0], p1[1], p1[2])
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p2[0], p2[1], p2[2])
    print >>f, 'gg::segEnd'
    print >>f, 'set %s [gg::conEnd]' % (label)
    return None

def rotateCopyConnector(f, label, p1, p2, angle, newLabel = None, replace = False):
    if replace is True:
        print >>f, 'gg::conTransformBegin $%s -maintain_linkage' % (label)
    else:
        print >>f, 'gg::conCopyBegin $%s' % (label)
    print >>f, 'gg::xformRotate {%f %f %f} {%f %f %f} %f' % (p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], angle)
    if replace is True:
        print >>f, 'gg::conTransformEnd'        
    else:
        print >>f, 'set %s [gg::conCopyEnd]' % (newLabel)
    return None

def rotateCopyDomain(f, label, p1, p2, angle, newLabel = None, replace = False):
    if replace is True:
        print >>f, 'gg::domTransformBegin $%s -maintain_linkage' % (label)
    else:
        print >>f, 'gg::domCopyBegin $%s' % (label)
    print >>f, 'gg::xformRotate {%f %f %f} {%f %f %f} %f' % (p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], angle)
    if replace is True:
        print >>f, 'gg::domTransformEnd'
    else:
        print >>f, 'set %s [gg::domCopyEnd]' % (newLabel)
    return None

def rotateCopyBlock(f, label, p1, p2, angle, newLabel = None, replace = False):
    if replace is True:
        print >>f, 'gg::blkTransformBegin $%s -maintain_linkage' % (label)
    else:
        print >>f, 'gg::blkCopyBegin $%s' % (label)
    print >>f, 'gg::xformRotate {%f %f %f} {%f %f %f} %f' % (p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], angle)
    if replace is True:
        print >>f, 'gg::blkTransformEnd'
    else:
        print >>f, 'set %s [gg::blkCopyEnd]' % (newLabel)        
    return None

def domain(f, edge1, edge2, edge3, edge4, label):
    print >>f, 'gg::domBegin -type STRUCTURED'
    print >>f, 'gg::edgeBegin'
    print >>f, 'gg::edgeAddCon $%s' % (edge1)
    print >>f, 'gg::edgeEnd'
    print >>f, 'gg::edgeBegin'
    print >>f, 'gg::edgeAddCon $%s' % (edge2)
    print >>f, 'gg::edgeEnd'
    print >>f, 'gg::edgeBegin'
    print >>f, 'gg::edgeAddCon $%s' % (edge3)
    print >>f, 'gg::edgeEnd'
    print >>f, 'gg::edgeBegin'
    print >>f, 'gg::edgeAddCon $%s' % (edge4)
    print >>f, 'gg::edgeEnd'
    print >>f, 'set %s [gg::domEnd]' % (label)
    return None

def writeGlyphFile(glyphFile, outputPrefix = 'OSUMach1.3', innerBlockSize = 0.5, innerBlockExponent = 1., 
                   nozzleRadius = 0.5, outerBlockRadius = 12.5, minAxialCoordinate = -10., maxAxialCoordinate = 34.,
                   potentialCoreLength = 6., nAxial = 1, nRadial = 281, nAzimuthal = 192, 
                   nInsideNozzleFraction = 0.1, nMinEquallySpaced = 12, minRadialSpacingFraction = 0.8, 
                   minAxialSpacingFraction = 2.25, nearInterfaceSpacingRatio = 1.):

    gridFile = 'Gridgen_' + outputPrefix + '_PLOT3D.x'

    assert nAzimuthal % 4 == 0

    nRadialInsideNozzle = math.floor(nInsideNozzleFraction * nRadial)
    radialSpacingNearInterface = nearInterfaceSpacingRatio * innerBlockSize / (nAzimuthal / 4)
    minRadialSpacing = minRadialSpacingFraction * innerBlockSize / (nAzimuthal / 4)
    minAxialSpacing = minAxialSpacingFraction * minRadialSpacing
    print minRadialSpacing

    f = open(glyphFile, 'w')

    print >>f, """
gg::memClear
gg::defReset
gg::tolReset
gg::dispViewReset
gg::dispGlide FALSE
set cwd [file dirname [info script]]
"""

    polarAngle = np.linspace(-np.pi / 4., np.pi / 4., 201)
    if innerBlockExponent < 2.:
        polarAngle += np.pi / 4.

    print >>f, 'gg::dbCurveBegin -type CUBIC'
    for theta in polarAngle:
        r = 0.5 * innerBlockSize / (abs(np.cos(theta)) ** innerBlockExponent + abs(np.sin(theta)) ** innerBlockExponent) ** (1. / innerBlockExponent)
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (r * np.cos(theta), r * np.sin(theta), 0.)
    print >>f, 'set db_inner_block_E [gg::dbCurveEnd]'
    print >>f, """
set inner_block_E [gg::conOnDBEnt $db_inner_block_E]
"""

    rotateCopyConnector(f, 'inner_block_E', [0., 0., 0.], [0., 0., 1.], 90., 'inner_block_N')
    rotateCopyConnector(f, 'inner_block_N', [0., 0., 0.], [0., 0., 1.], 90., 'inner_block_W')
    rotateCopyConnector(f, 'inner_block_W', [0., 0., 0.], [0., 0., 1.], 90., 'inner_block_S')

    # Dimension inner block and create domain.
    print >>f, 'gg::conDim [list $inner_block_E $inner_block_N $inner_block_W $inner_block_S] %i' % (nAzimuthal / 4 + 1)
    domain(f, 'inner_block_S', 'inner_block_E', 'inner_block_N', 'inner_block_W', 'dom_inner')
    if innerBlockExponent < 2.:
        rotateCopyDomain(f, 'dom_inner', [0., 0., 0.], [0., 0., 1.], -45., replace = True)

    # Create south-east connector to specify extrusion step size.
    polarAngleSE = polarAngle[0]
    radialCoordinateSE = 0.5 * innerBlockSize / (abs(np.cos(polarAngleSE)) ** innerBlockExponent + abs(np.sin(polarAngleSE)) ** innerBlockExponent) ** (1. / innerBlockExponent)

    drawLine3D(f, 'tmp_extrusion_SE',
               [radialCoordinateSE * np.cos(polarAngleSE), radialCoordinateSE * np.sin(polarAngleSE), 0.],
               [outerBlockRadius * np.cos(polarAngleSE), outerBlockRadius * np.sin(polarAngleSE), 0.])

    # Dimension south-east connector.
    print >>f, 'gg::conDim [list $tmp_extrusion_SE] %i' % (nRadial)
    print >>f, 'gg::conSetBreakPt $tmp_extrusion_SE {%f %f 0}' % ((radialCoordinateSE + nMinEquallySpaced * radialSpacingNearInterface) * np.cos(polarAngleSE), (radialCoordinateSE + nMinEquallySpaced * radialSpacingNearInterface) * np.sin(polarAngleSE))
    print >>f, 'gg::conSetBreakPt $tmp_extrusion_SE {%f %f 0}' % (nozzleRadius * np.cos(polarAngleSE), nozzleRadius * np.sin(polarAngleSE)) # radial location of nozzle
    print >>f, 'gg::conSubConDim $tmp_extrusion_SE [list %i %i %i]' % (nMinEquallySpaced + 1, nRadialInsideNozzle - nMinEquallySpaced + 1, nRadial - nRadialInsideNozzle)
    print >>f, 'gg::conBeginSpacing $tmp_extrusion_SE -sub 3 %f' % (minRadialSpacing)
    print >>f, 'gg::conBeginSpacing $tmp_extrusion_SE -sub 2 %f' % (radialSpacingNearInterface)
    print >>f, 'gg::conEndSpacing $tmp_extrusion_SE -sub 2 %f' % (minRadialSpacing)
    print >>f, 'gg::conBeginSpacing $tmp_extrusion_SE -sub 1 %f' % (radialSpacingNearInterface)
    print >>f, 'gg::conEndSpacing $tmp_extrusion_SE -sub 1 %f' % (radialSpacingNearInterface)
    if innerBlockExponent < 2.:
        rotateCopyConnector(f, 'tmp_extrusion_SE', [0., 0., 0.], [0., 0., 1.], -45., replace = True)

    # Extrude to create the outer blocks.
    print >>f, """
gg::domExtrusionBegin [list $inner_block_E $inner_block_N $inner_block_W $inner_block_S] -default HYPERBOLIC
gg::domExtrusionAtt -flip
gg::domExtrusionAtt -growth_subcon [list [list $tmp_extrusion_SE 1] [list $tmp_extrusion_SE 2] [list $tmp_extrusion_SE 3]]
"""
    print >>f, 'gg::domExtrusionAtt -stop_height %f' % (outerBlockRadius)
    print >>f, 'gg::domExtrusionStep %i' % (nRadial)
    print >>f, """
set dom_outer_blocks [gg::domExtrusionEnd]
set dom_outer_E [lindex $dom_outer_blocks 0]
set dom_outer_N [lindex $dom_outer_blocks 1]
set dom_outer_W [lindex $dom_outer_blocks 2]
set dom_outer_S [lindex $dom_outer_blocks 3]
unset dom_outer_blocks
"""

    if nAxial > 1:

        # Move domains to axial start coordinate.
        print >>f, 'gg::domTransformBegin [list $dom_inner $dom_outer_E $dom_outer_N $dom_outer_W $dom_outer_S] -maintain_linkage'
        print >>f, 'gg::xformTranslate {0 0 %f}' % (minAxialCoordinate)
        print >>f, 'gg::domTransformEnd'

        # Create the axis connector.
        drawLine3D(f, 'axis',
                   [0., 0., minAxialCoordinate],
                   [0., 0., maxAxialCoordinate])
        print >>f, 'gg::conDim [list $axis] %i' % (nAxial)
        print >>f, 'gg::conSetBreakPt $axis {0 0 %f}' % (potentialCoreLength)
        print >>f, 'gg::conBreakPtSpacing $axis -breakpt 1 %f' % (minAxialSpacing)

        # Extrude blocks.
        print >>f, 'gg::blkExtrusionBegin [list $dom_inner $dom_outer_E $dom_outer_N $dom_outer_W $dom_outer_S] -default PATH'
        print >>f, 'gg::blkExtrusionAtt -path_connector [list [list $axis 1] [list $axis 2]]'
        print >>f, 'gg::blkExtrusionStep %i' % (nAxial)
        print >>f, """
set all_blocks [gg::blkExtrusionEnd]
set blk_inner [lindex $all_blocks 0]
set blk_outer_E [lindex $all_blocks 1]
set blk_outer_N [lindex $all_blocks 2]
set blk_outer_W [lindex $all_blocks 3]
set blk_outer_S [lindex $all_blocks 4]
unset all_blocks
"""

        # Re-orient blocks.
        print >>f, 'gg::blkSpecifyIJK [list $blk_outer_E $blk_outer_N $blk_outer_W $blk_outer_S] 5 4 1'

        # Save PLOT3D grid file.
        print >>f, 'gg::blkExport ALL [file join $cwd "%s"] -style PLOT3D -form UNFORMATTED -precision DOUBLE -endian NATIVE' % (gridFile)

    else:

        # Save PLOT3D grid file.
        print >>f, 'gg::domExport ALL [file join $cwd "%s"] -style PLOT3D -form UNFORMATTED -precision DOUBLE -endian NATIVE' % (gridFile)

    f.close()

if __name__ == '__main__':

    glyphFile = 'gridgenScript.glf'
    writeGlyphFile(glyphFile, innerBlockExponent = 1.12148531779, nAxial = 481, potentialCoreLength = 11., nearInterfaceSpacingRatio = 0.999793250757)
    subprocess.check_output(["gridgen", "-b", glyphFile], stderr = subprocess.STDOUT)
