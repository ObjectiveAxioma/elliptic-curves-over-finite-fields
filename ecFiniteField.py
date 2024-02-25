import numpy as np

##########################################
### ELLIPTIC CURVES OVER FINITE FIELDS ###
##########################################

##########################################################################################
# This program is motivated by the fact that my previous EC implementation cannot handle #
# finite fields conveniently and I'd like to implement some algorithms I'm learning for  #
# them to improve my understanding and learn to work with finite fields.                 # 
##########################################################################################



##########################################################################################
# The Point class holds a point on the curve.                                            #
#                                                                                        #
# Point.x: the x-coordinate of the point in Fq                                           #
# Point.y: the y-coordinate of the point in Fq                                           #
# Point.curve: the curve over which the point is defined                                 #
#                                                                                        #
# Point addition, negation, subtraction, and mult by integer are supported with magic    #
# methods. Operations will be performed mod q when necessary.                            #
##########################################################################################
class Point():

    def __init__(self, x, y, E):
        if isinstance(E, EllipticCurve):
            self.curve = E
        else:
            return "A curve must be provided with P."
        self.x = x
        self.y = y

    def __neg__(self):
        if self.x == float("inf") or self.y == float("inf"):
            return Point(self.x, self.y, self.curve)
        return Point(self.x, self.curve.char - self.y, self.curve)

    def __add__(self, Q):
        return self.curve.addPoints(self, Q)

    def __sub__(self, Q):
        return self.curve.addPoints(self, -Q)

    def __mul__(self, n):
        if isinstance(n, int):
            return self.curve.mulPoint(self, n)

    __lmul__ = __mul__
    __rmul__ = __mul__

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __str__(self):
        return str("(" + str(self.x) + ", " + str(self.y) + ")")

    __repr__ = __str__



##########################################################################################
# The EllipticCurve class represents the elliptic curve over the finite field Fq.        #
#                                                                                        #
# EllipticCurve.char: characteristic of Fq over which the curve is defined               #
# EllipticCurve.q: order of the finite field Fq over which the curve is defined          #
# EllipticCurve.a: first coefficient of Weierstrass or general Weierstrass if char = 2   #
# EllipticCurve.b: second coefficient of Weierstrass or  general Weierstrass is char = 2 #
# EllipticCurve.c: the third coefficient of gen. Weierst., only when char = 2            #
# EllipticCurve.dis: discriminant of the curve, must be != 0                             #
#                                                                                        #
# isOnCurve(P): returns 1 if P is on curve, 0 otherwise                                  #
# addPoints(P, Q): adds points mod char                                                  #
# timesTwo(P): doubles the point P mod char                                              #
# mulPoint(n, P): adds P to itself n times mod char                                      #
# plotCurve(): plots the curve on the affine 2-plane of Fq                               #
#                                                                                        #
# Comparison with "==" is supported. Calling the object returns (general) Weierst. eqn   #
##########################################################################################
class EllipticCurve():

    # Initializes an EllipticCurve object with characteristic and coefficients
    def __init__(self, char, q, a, b, c=None):
        self.char = char
        self.q = q
        self.a = a
        self.b = b
        if self.char == 2 and (isinstance(c, int) or isinstance(c, float)):
            self.c = c

        self.eqn = "w"
        if self.char == 2:
            # need to add cases for a0 = 0 and a0 != 0


