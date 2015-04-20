#!/usr/bin/env python

def derivative(a, order = 1, axis = -1, scheme = 'SBP 1-2', periodic = False):

    import numpy as np
    from math import copysign

    """
    Calculate the n-th order finite-difference along given axis using a diagonal-norm summation-by-parts (SBP) operator.

    Parameters
    ----------
    a: array_like
       Input array
    order: {1, 2}, optional
           The order of the SBP operator, default is 1.
    axis: int, optional
          The axis along which the finite-difference is calculated, default is the last axis.
    scheme: {'SBP 1-2', 'SBP 2-4', 'SBP 3-6', 'SBP 4-8'}, optional
            The SBP scheme used to calculate the finite-difference, default is 'SBP 1-2'.

    Returns
    -------
    derivative: ndarray
                An array of the same shape as `a` giving the finite-difference approximation of `a` using the SBP operator given by `scheme`.

    """
    
    allowedOrders = [1, 2]
    if order not in allowedOrders:
        raise ValueError("'%s' is not a valid order (allowed values are: %s)" % (order, ", ".join(map(str, allowedOrders))))

    allowedSchemes = ['SBP 1-2', 'SBP 2-4', 'SBP 3-6', 'SBP 4-8']
    if scheme not in allowedSchemes:
        raise ValueError("'%s' is not a valid scheme (allowed values are: %s)" % (scheme, ", ".join(map(str, allowedSchemes))))

    a = np.asanyarray(a)
    b = np.empty_like(a)

    if order == 1:
        if scheme == 'SBP 1-2':
            rhsInterior = np.array([0., 1./2.])
            rhsBoundary = np.zeros([2, 1])
            rhsBoundary[:,0] = [-1., 1.]
        elif scheme == 'SBP 2-4':
            rhsInterior = np.array([0., 2./3., -1./12.])
            rhsBoundary = np.zeros([6, 4])
            rhsBoundary[:4,0] = [-24./17., 59./34., -4./17., -3./34.]
            rhsBoundary[:3,1] = [-1./2., 0., 1./2.]
            rhsBoundary[:5,2] = [4./43., -59./86., 0., 59./86., -4./43.]
            rhsBoundary[:6,3] = [3./98., 0., -59./98., 0., 32./49., -4./49.]
        elif scheme == 'SBP 3-6':
            rhsInterior = np.array([0., 3./4., -3./20., 1./60.])
            rhsBoundary = np.zeros([9, 6])
            rhsBoundary[:6,0] = [-21600./13649., 104009./54596., 30443./81894., -33311./27298., 16863./27298., -15025./163788.]
            rhsBoundary[:6,1] = [-104009./240260., 0., -311./72078., 20229./24026., -24337./48052., 36661./360390.]
            rhsBoundary[:6,2] = [-30443./162660., 311./32532., 0., -11155./16266., 41287./32532., -21999./54220.]
            rhsBoundary[:7,3] = [33311./107180., -20229./21436., 485./1398., 0., 4147./21436., 25427./321540., 72./5359.]
            rhsBoundary[:8,4] = [-16863./78770., 24337./31508., -41287./47262., -4147./15754., 0., 342523./472620., -1296./7877., 144./7877.]
            rhsBoundary[:9,5] = [15025./525612., -36661./262806., 21999./87602., -25427./262806., -342523./525612., 0., 32400./43801., -6480./43801., 720./43801.]
        else:
            rhsInterior = np.array([0., 4./5., -1./5., 4./105., -1./280.])
            x1 = 541./1000.
            x2 = -27./400.
            x3 = 187./250.
            rhsBoundary = np.zeros([12, 8])
            rhsBoundary[:8,0] = [-2540160./1498139., 9.*(2257920.*x1+11289600.*x2+22579200.*x3-15849163.)/5992556., 3.*(-33868800.*x1-162570240.*x2-304819200.*x3+235236677.)/5992556., (609638400.*x1+2743372800.*x2+4572288000.*x3-3577778591.)/17977668., 3.*(-16934400*x1-67737600.*x2-84672000.*x3+67906303.)/1498139., 105.*(967680.*x1+2903040.*x2-305821.)/5992556., 49.*(-1244160.*x1+18662400.*x3-13322233.)/17977668., 3.*(-6773760.*x2-33868800.*x3+24839327.)/5992556.]
            rhsBoundary[:8,1] = [9.*(-2257920.*x1-11289600.*x2-22579200.*x3+15849163.)/31004596., 0., 3.*(7257600.*x1+33868800.*x2+60963840.*x3-47167457.)/2214614., 3.*(-9676800.*x1-42336000.*x2-67737600.*x3+53224573.)/1107307., 7.*(55987200.*x1+217728000.*x2+261273600.*x3-211102099.)/13287684., 3.*(-11612160.*x1-33868800.*x2+3884117.)/2214614., 150.*(24192.*x1-338688.*x3+240463.)/1107307., (152409600.*x2+731566080.*x3-536324953.)/46506894.]
            rhsBoundary[:8,2] = [(33868800.*x1+162570240.*x2+304819200.*x3-235236677.)/1743924., (-7257600.*x1-33868800.*x2-60963840.*x3+47167457.)/124566., 0., (24192000.*x1+101606400.*x2+152409600.*x3-120219461.)/124566., (-72576000.*x1-270950400.*x2-304819200.*x3+249289259.)/249132., 9.*(806400.*x1+2257920.*x2-290167.)/41522., 6.*(-134400.*x1+1693440.*x3-1191611.)/20761., 5.*(-2257920.*x2-10160640.*x3+7439833.)/290654.]
            rhsBoundary[:8,3] = [(-609638400.*x1-2743372800.*x2-4572288000.*x3+3577778591.)/109619916., 3.*(9676800.*x1+42336000.*x2+67737600.*x3-53224573.)/1304999., 3.*(-24192000.*x1-101606400.*x2-152409600.*x3+120219461.)/2609998., 0., 9.*(16128000.*x1+56448000.*x2+56448000.*x3-47206049.)/5219996., 3.*(-19353600.*x1-50803200.*x2+7628371.)/2609998., 2.*(10886400.*x1-114307200.*x3+79048289.)/3914997., 75.*(1354752.*x2+5419008.*x3-3952831.)/18269986.]
            rhsBoundary[:9,4] = [3.*(16934400.*x1+67737600.*x2+84672000.*x3-67906303.)/2096689., 7.*(-55987200.*x1-217728000.*x2-261273600.*x3+211102099.)/3594324., 3.*(72576000.*x1+270950400.*x2+304819200.*x3-249289259.)/1198108., 9.*(-16128000.*x1-56448000.*x2-56448000.*x3+47206049.)/1198108., 0., 105.*(414720.*x1+967680.*x2-165527.)/1198108., 15.*(-967680.*x1+6773760.*x3-4472029.)/1198108., (-304819200.*x2-914457600.*x3+657798011.)/25160268., -2592./299527.]
            rhsBoundary[:10,5] = [5.*(-967680.*x1-2903040.*x2+305821.)/1237164., (11612160.*x1+33868800.*x2-3884117.)/618582., 9.*(-806400.*x1-2257920.*x2+290167.)/206194., (19353600.*x1+50803200.*x2-7628371.)/618582., 35.*(-414720.*x1-967680.*x2+165527.)/1237164., 0., 80640.*x1/103097., 80640.*x2/103097., 3072./103097., -288./103097.]
            rhsBoundary[:11,6] = [7.*(1244160.*x1-18662400.*x3+13322233.)/8041092., 150.*(-24192.*x1+338688.*x3-240463.)/670091., 54.*(134400.*x1-1693440.*x3+1191611.)/670091., 2.*(-10886400.*x1+114307200.*x3-79048289.)/2010273., 15.*(967680.*x1-6773760.*x3+4472029.)/2680364., -725760.*x1/670091., 0., 725760.*x3/670091., -145152./670091., 27648./670091., -2592./670091.]
            rhsBoundary[:12,7] = [3.*(6773760.*x2+33868800.*x3-24839327.)/20510956., (-152409600.*x2-731566080.*x3+536324953.)/30766434., 45.*(2257920.*x2+10160640.*x3-7439833.)/10255478., 75.*(-1354752.*x2-5419008.*x3+3952831.)/10255478., (304819200.*x2+914457600.*x3-657798011.)/61532868., -5080320.*x2/5127739., -5080320.*x3/5127739., 0., 4064256./5127739., -1016064./5127739., 193536./5127739., -18144./5127739.]

    nd = len(a.shape)
    n = a.shape[axis]

    width = rhsBoundary.shape[0]
    depth = rhsBoundary.shape[1]
    interiorSize = rhsInterior.size

    slice1 = [slice(None)] * nd
    slice2 = [slice(None)] * nd
    slice3 = [slice(None)] * nd
    slice1[axis] = slice(depth, -depth)

    b[slice1] = 0.
    if order % 2 == 0:
        for i in range(interiorSize):
            slice2[axis] = slice(depth + i, n - depth + i)
            slice3[axis] = slice(depth - i, n - depth - i)
            b[slice1] += rhsInterior[i] * (a[slice2] + a[slice3])
    else:
        for i in range(1, interiorSize):
            slice2[axis] = slice(depth + i, n - depth + i)
            slice3[axis] = slice(depth - i, n - depth - i)
            b[slice1] += rhsInterior[i] * (a[slice2] - a[slice3])

    slice1[axis] = slice(None, depth)
    slice2[axis] = slice(None, width)
    b[slice1] = np.rollaxis(np.tensordot(rhsBoundary, a[slice2], axes = (0, axis)), axis)

    slice1[axis] = slice(None, -depth - 1, -1)
    slice2[axis] = slice(None, -width - 1, -1)
    if order % 2 == 0:
        b[slice1] = np.tensordot(a[slice2], rhsBoundary, axes = (axis, 0))
    else:
        b[slice1] = - np.rollaxis(np.tensordot(rhsBoundary, a[slice2], axes = (0, axis)), axis)

    return b
