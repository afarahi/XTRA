from math import sin, cos, tan, atan, asin, acos, atan2
from XCat_Objects import pi, DtoR, D2R, R2D
import numpy as np
from numpy.linalg import norm


# RA  = phi
# DEC = pi/2 - theta (theta = pi/2 - DEC)
def DEC2theta(DEC):
    return np.pi/2. - DEC


def theta2DEC(theta):
    return np.pi/2. - theta


def Cartesian2Spherical(x,y,z):
    cart_vector = [x, y, z]
    r = norm(cart_vector)
    unit = cart_vector / r
    theta = acos(unit[2])
    phi = atan2(unit[1], unit[0])
    return r, theta, phi


def Spherical2Cartesian(r,RA,DEC):
    cos_DEC = cos(DEC)
    x = r * cos(RA) * cos_DEC
    y = r * sin(RA) * cos_DEC
    z = r * sin(DEC)
    return x, y, z


def Cartesian2RADEC(x,y,z,noNeg=True):
    cart_vector = [x, y, z]
    r = norm(cart_vector)  # it is not checking for r=1.0 but it should
    unit = cart_vector / r
    theta = acos(unit[2])
    phi = atan2(unit[1], unit[0])  # RA = phi, DEC = theta2DEC(theta)
    if noNeg:
        while (phi < 0.0): phi += 2. * np.pi
    return phi, theta2DEC(theta)


# RA_o and DEC_o such that RA_o and DEC_o 
# in spherical coordinate become zeros 
def CoordTrans(x,y,z,RA_o,DEC_o):
    cosD = cos(DEC_o); sinD = sin(DEC_o)
    cosR = cos(RA_o);  sinR = sin(RA_o)
    xn =  cosD*cosR*x + cosD*sinR*y + sinD*z
    yn =      -sinR*x +      cosR*y
    zn = -sinD*cosR*x - sinD*sinR*y + cosD*z
    return xn,yn,zn


def invCoordTrans(x,y,z,RA_o,DEC_o):
    cosD = cos(DEC_o); sinD = sin(DEC_o)
    cosR = cos(RA_o);  sinR = sin(RA_o)
    xn =  cosD*cosR*x - sinR*y - cosR*sinD*z
    yn =  sinR*cosD*x + cosR*y - sinR*sinD*z
    zn =    sinD*x             +   cosD*z
    return xn,yn,zn


# RA, DEC befor rotation (in degree)
# RA_o, DEC_o the new origin for coordinate (in degree)
# return new RA, DEC (in degree)
def rotatingCoord(RA,DEC,RA_o,DEC_o,noNeg=True):
    x,y,z = Spherical2Cartesian(1.0,D2R*RA,D2R*DEC)
    x,y,z = CoordTrans(x,y,z,D2R*RA_o,D2R*DEC_o)
    RAn,DECn = Cartesian2RADEC(x,y,z,noNeg=noNeg)
    return R2D*RAn,R2D*DECn


def EquatorialTOGalactic(RA,DEC):
    # l
    gLong = 303.0*DtoR - atan( sin(192.25*DtoR-RA*DtoR) / (cos(192.25*DtoR-RA*DtoR)*sin(27.4*DtoR) - tan(DEC*DtoR)*cos(27.4*DtoR)) )
    # b
    gLati = asin(sin(DEC*DtoR)*sin(27.4*DtoR) + cos(DEC*DtoR)*cos(27.4*DtoR)*cos(192.25*DtoR-RA*DtoR))
    return (gLong/DtoR,gLati/DtoR)


''' EXAMPLES
#EXAMPLE (WORKING)
#print "RA & DEC :", rotatingCoord(20.0,30.0,20.0,30.0)
#RA_o  = 00.*D2R
#DEC_o = 30.*D2R
#RA  = 00.*D2R
#DEC = 20.*D2R
#x,y,z = Spherical2Cartesian(1.0,RA,DEC)
#x,y,z = CoordTrans(x,y,z,RA_o,DEC_o)
#RA,DEC = Cartesian2RADEC(x,y,z)
#print "RA  : ", R2D*RA; print "DEC : ", R2D*DEC

#x=4.0; y=1.0; z=3.0
#x,y,z = CoordTrans(x,y,z,1.0,1.0)
#x,y,z = invCoordTrans(x,y,z,1.0,1.0)
#print x,y,z


RA_o  = 30.*D2R
DEC_o = 20.*D2R
r     = 1.0

RA  = 50.*D2R
DEC = 20.*D2R
r   = 1.0

print 'RA  = ', R2D*RA
print 'DEC = ', R2D*DEC
print 'r   = ', r

x,y,z = Spherical2Cartesian(r,RA,DEC)
x,y,z = CoordTrans(x,y,z,RA_o,DEC_o)
r, theta, phi = Cartesian2Spherical(x,y,z)

print '--------------------------------'
print 'RA  = ', R2D*phi
print 'DEC = ', R2D*theta2DEC(theta)
print 'r   = ', r
'''