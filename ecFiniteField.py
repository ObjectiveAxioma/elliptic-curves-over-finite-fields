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
        self.x = x % E.char
        self.y = y % E.char

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
# EllipticCurve.c: third coefficient of gen. Weierst., only when char = 3 or 2 and t = 1 #
# EllipticCurve.t: whether or not a1 = 0 in gen. Weierst., only when char = 2            #
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
    def __init__(self, char, q, a, b, c=None, t=None):
        self.char = char
        self.q = q
        self.a = a % char
        self.b = b % char
        if isinstance(t, int):
            if t == 1 or t == 0:
                self.t = t
        if self.char == 2 and (isinstance(c, int) or isinstance(c, float)):
            self.eqn = "n"          # n for non-zero, i.e. a1 != 0
            if self.t == 0:
                self.c = c % 2 
                self.eqn = "z"      # z for zero, i.e. a1 = 0

        if self.char != 2:
            self.eqn = "w"

        if self.char == 3:
            self.eqn = "t"      # t for three, i.e. char = 3
            self.c = c % 3

    def isOnCurve(self, P):
        if self.char != 2 and self.char != 3:
            toCheck = (P.y)**2 - (P.x)**3 - self.a*(P.x) - self.b
            if toCheck == 0:
                return 1
            else:
                return 0
        elif self.char == 3:
            toCheck = (P.y)**2 - (P.x)**3 - self.a*(P.x)**2 - self.b*(P.x) - self.c

    # Finds the order of E over the field Fp using Legendre symbols and then
    # uses the fact that q = p^k implies |E(Fq)| = q^k + 1 - (alpha + beta)^k
    def orderOfCurve(self):
        # Get the Legendre symbols over Fp
        symbols = getLegendre(self.char)
        print("symbols = " + str(symbols) + "\n")
        primeOrder = self.char + 1                      # The order of Fp to be used to find |Fq|

        # Loop to find the order of Fp from the Legendre symbols
        n = 0
        while n < self.char:
            primeOrder += symbols[int(n*n*n + self.a*n + self.b) % self.char]
            n+=1

        # Finds k for which q = p^k
        k = int(np.log(self.q) / np.log(self.char))
        print("k = " + str(k) + "\n")

        # Generate the values of s_n = (alpha + beta)^n with alpha and beta being roots of the 
        # characteristic polynomial of phi_n; necessary for finding the order of Fq
        x0 = 2
        x1 = self.char + 1 - primeOrder
        xvals = [x0, x1]
        n = 2
        while n < k:
            xvals.append(x1 * xvals[n] - self.char * xvals[n-1])
            n+=1
        
        x = xvals[len(xvals) - 1]
        if k == 1:
            x = x0
        print("x = " + str(x) + "\n")
        if k > 1:
            order = (self.q)^k + 1 - x^k
        else:
            order = primeOrder

        return order

    def orderOfPoint(self, P):
        return 0

    def __str__(self):
        # When characteristic is 2
        if self.char == 2:
            # a1 != 0
            if self.eqn == 'n':
                return str("y^2 + xy = x^3 + " + str(self.a) + "x^2 + " + str(self.b) + ", char = " + str(self.char) + ", field = F" + str(self.q))
            
            # a1 = 0
            if self.eqn == 'z':
                return str("y^2 + " + str(self.a) + "y = x^3 + " + str(self.b) + "x + " + str(self.c) + ", char = " + str(self.char) + ", field = F" + str(self.q))

        # When characteristic is 3
        if self.char == 3:
            return str("y^2 = x^3 + " + str(self.a) + "x^2 + " + str(self.b) + "x + " + str(self.a) + ", char = " + str(self.char) + ", field = F" + str(self.q))

        # When characteristic is not 2 or 3
        if self.a > 0 and self.b > 0:
            return str("y^2 = x^3 + " + str(self.a) + "x + " + str(self.b) + ", char = " + str(self.char) + ", field = F" + str(self.q))
        elif self.a > 0 and self.b < 0:
            return str("y^2 = x^3 + " + str(self.a) + "x - " + str(-self.b) + ", char = " + str(self.char) + ", field = F" + str(self.q))
        elif self.a < 0 and self.b > 0:
            return str("y^2 = x^3 - " + str(-self.a) + "x + " + str(self.b) + ", char = " + str(self.char) + ", field = F" + str(self.q))
        elif self.a < 0 and self.b < 0:
            return str("y^2 = x^3 - " + str(-self.a) + "x - " + str(-self.b) + ", char = " + str(self.char) + ", field = F" + str(self.q))

    __repr__ = __str__


# Helper method to compute the Ledengre symbold of all numbers from 1 to p-1. This is used
# to find the order of the curve over Fq by first finding the order of the curve over Fp.
# where q = p^n, and then applying |E(q^n)| = q^n + 1 - (alpha^n + beta^n) for roots of the
# characteristic equation of the Frobenius endomorphism.
# Note: assume char != 2, 3.
def getLegendre(p):
    
    # Get a list of all squares and use it to compute the Legendre symbols of all elements in Fq.
    n = 2
    squares = []
    while n < p:
        if n not in squares:
            squares.append(n*n)
        n+=1

    # Get a list of all the legendre symbols of elemnts of Fq.
    n = 2
    legendre = [0, 1]
    while n < p:
        if n in squares:
            legendre.append(1)
        elif n not in squares:
            legendre.append(-1)
        n+=1

    return legendre




#### Current testing code ####
curve = EllipticCurve(5, 5, 1, 1)       # y^2 = x^3 + x + 1 over F5
order = curve.orderOfCurve()
curve2 = EllipticCurve(5, 25, 1, 1)     # y^2 = x^3 + x + 1 over F25
order2 = curve2.orderOfCurve()

print("Order of y^2 = x^3 + x + 1 over F5: ")
print(order)
print("\n")
print("Orer of y^2 - x^3 _ x _ 1 over F25: ")
print(order2)
