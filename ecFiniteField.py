import numpy as np

##########################################
### ELLIPTIC CURVES OVER FINITE FIELDS ###
##########################################


##########################################################################################
# This program contains the Point and EllipticCurve classes, constituting a full         #
# implementation of the elliptic curve object over a finite field with algorithms for    #
# point addition, multiplication with integers, finding group and point order, and more. #
# This is to be used in my future projects in elliptic curve cryptography.               #
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
        if x == float("inf"):
            self.x = float("inf")
        if y == float("inf"):
            self.y = float("inf")

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
        if other == "inf":
            return self.x == float("inf") and self.y == float("inf")
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
# EllipticCurve.j: the j-invariant of the curve                                          #
# EllipticCurve.dis: discriminant of the curve, must be != 0                             #
#                                                                                        #
# isOnCurve(P): returns 1 if P is on curve, 0 otherwise                                  #
# isIsomorphic(E) : returns 1 if E and C are isomorphic and 0 otherwise                  #
# timesTwo(P): doubles the point P mod char                                              #
# addPoints(P, Q): adds P and Q according to the group law on E mod char of Fq           #
# mulPoint(n, P): adds P to itself n times mod char                                      #
# orderOfCurve(): computes the order of the curve over Fq using Lagrange symbols         #
# orderOfPoint(P): computes the order of a point using baby step, giant step             #
# isSupersingular(): returns 0 if E[p] != 0, otherwise returns 1                         #
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

        self.j = (4*(self.a)**3) / (4*(self.a)**3 + 27*(self.b)**2)

    # Check if a given  point is on the curve.
    def isOnCurve(self, P):
        # The point at infinity is on every curve.
        if P.x == float("inf") or P.y == float("inf"):
            return 1;

        # If the character of the curve is not 2 or 3, follow standard procedure.
        if self.char != 2 and self.char != 3:
            toCheck = ((P.y)**2 - (P.x)**3 - self.a*(P.x) - self.b) % self.char
            if toCheck == 0:
                return 1
            else:
                return 0

        # If the character is 3, another formula must be used.
        elif self.char == 3:
            toCheck = (P.y)**2 - (P.x)**3 - self.a*(P.x)**2 - self.b*(P.x) - self.c % 3
            if toCheck == 0:
                return 1
            else:
                return 0

    # Determines whether the curve is isomorphic to another curve E by determining whether
    # or not their j-invariants are identical.
    def isIsomorphic(self, E):
        if self.j == E.j:
            return 1

        return 0

    # Adds a point of the curve to itself over the finite field of definition.
    def timesTwo(self, P):
        # Find the inverses of all values in Fp.
        inverses = findInverses(self.char)

        m = ((3 * P.x * P.x + self.a) * inverses[2] * inverses[P.y % self.char]) % self.char

        x = (m*m - 2*P.x) % self.char
        y = -(m*(x - P.x) + P.y)
        Q = Point(x, y, self)

        return Q

    # Add the points P and Q according to the group law on the curve E, modding the cooridnates
    # by the characteristic of the curve.
    def addPoints(self, P, Q):
        if not isinstance(P, Point) or not isinstance(Q, Point):
            return "Please enter a Point object."
        if not self.isOnCurve(P) or not self.isOnCurve(Q):
            return "Please enter a Point on the curve."

        # Find all of the inverses mod the char of the elliptic curve so that the point
        # addition formula can be used.
        inverses = findInverses(self.char)

        if not P.x == float("inf") and not P.y == float("inf") and not Q.x == float("inf") and not Q.y == float("inf"):
            m = (Q.y - P.y) * inverses[int(Q.x - P.x) % self.char]
            x = m*m - Q.x - P.x % self.char
            y = -(m*(x - P.x) + P.y) % self.char

        if P.x == Q.x:
            if P.y != Q.y:
                return Point(float("inf"), float("inf"), self)
            else:
                return self.timesTwo(P)

        if P.x == float("inf") and P.y == float("inf"):
            x = Q.x
            y = Q.y

        if Q.x == float("inf") and Q.y == float("inf"):
            x = P.x
            y = P.y

        R = Point(x, y, self)
        return R

    # Multiply a point by an integer n, i.e. add it to itself n times.
    def mulPoint(self, P, n):
        if P.x == float("inf") or P.y == float("inf"):
            return P

        # Deal with the case where n is negative.
        sign = 1
        if n < 0:
            n = -n
            sign = -1


        a = n
        B = Point(float("inf"), float("inf"), self)
        C = P

        while a != 0:
            if a %2 == 0:
                a = a/2
                B = B
                C = self.timesTwo(C)
            elif a % 2 == 1:
                a = a - 1
                B = B + C
                C = C
    
        if sign == -1:
            return -B
        return B

    # Finds the order of E over the field Fp using Legendre symbols and then
    # uses the fact that q = p^k implies |E(Fq)| = q^k + 1 - (alpha + beta)^k
    def orderOfCurve(self):
        # Get the Legendre symbols over Fp
        symbols = getLegendre(self.char)
        primeOrder = self.char + 1                      # The order of Fp to be used to find |Fq|

        # Loop to find the order of Fp from the Legendre symbols
        n = 0
        while n < self.char:
            primeOrder += symbols[int(n*n*n + self.a*n + self.b) % self.char]
            n+=1

        # Finds k for which q = p^k
        k = int(np.log(self.q) / np.log(self.char))

        # Generate the values of s_n = (alpha + beta)^n with alpha and beta being roots of the 
        # characteristic polynomial of phi_n; necessary for finding the order of Fq
        x0 = 2
        x1 = self.char + 1 - primeOrder
        xvals = [x0, x1]
        n = 1
        while n < k:
            xvals.append(x1 * xvals[n] - self.char * xvals[n-1])
            n+=1
        
        x = xvals[len(xvals) - 1]
        if k == 1:
            x = x1
        if k > 1:
            order = (self.char)**k + 1 - x
        else:
            order = primeOrder

        return order

    # Finds the order of a point using the baby step, giant step algorithm.
    def orderOfPoint(self, P):
        Q = (self.q + 1) * P
        m = int(self.q**(1/4) + 1)

        k = -m
        isMatch = 0
        match = 0
        sign = 0
        while k <= m and isMatch == 0:
            x = Q + k*(2*m*P)
            j = 0
            while j <+ m and isMatch == 0:
                currentPoint = j*P
                if x == currentPoint:
                    isMatch = 1
                    sign = -1
                elif x == -currentPoint:
                    isMatch = 1
                    sign = 1
                j+=1
            k+=1

        M = self.q + 1 + 2*m*(k-1) + (sign) * (j-1)
        M = int(((M**2)**(1/2)))

        return self.recursiveBSGS(M, P)

    # Computes the recursive part of the baby step, giant step algorithm.
    def recursiveBSGS(self, M, P):
        print("M = " + str(M))
        factors = factor(M)
        repeat = 0
        x = 0
        j = 1
        while j < len(factors):
            x = int(M/factors[j])
            if x*P == "inf":
                repeat = x
            j+=1
        if repeat != 0:
            return self.recursiveBSGS(repeat, P)
        return M


    # Determines whether or not the curve is supersingular using the equivalent condition
    # that E/Fp is supersingular if and only if the order of the curve is 1 mod char.
    def isSupersingular(self):
        order = self.orderOfCurve()
        if order % self.char == 1:
            return 1
        
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
            if self.a > 0 and self.b > 0 and self.c > 0:
                return str("y^2 = x^3 + " + str(self.a) + "x^2 + " + str(self.b) + "x + " + str(self.c) + ", char = " + str(self.char) + ", field = F" + str(self.q))
            elif self.a > 0 and self.b > 0 and self.c == 0:
                return str("y^2 = x^3 + " + str(self.a) + "x^2 + " + str(self.b) + "x" + ", char = " + str(self.char) + ", field = F" + str(self.q))
            elif self.a > 0 and self.b == 0 and self.c > 0:
                return str("y^2 = x^3 + " + str(self.a) + "x^2 + " + str(self.c) +", char = " + str(self.char) + ", field = F" + str(self.q))
            elif self.a == 0 and self.b > 0 and self.c > 0:
                return str("y^2 = x^3 + " + str(self.b) + "x + " + str(self.c) + ", char = " + str(self.char) + ", field = F" + str(self.q))
            elif self.a > 0 and self.b == 0 and self.c == 0:
                return str("y^2 = x^3 + " + str(self.a) + ", char = " + str(self.char) + ", field = F" + str(self.q))
            elif self.a == 0 and self.b == 0 and self.c > 0:
                return str("y^2 = x^3 + " + str(self.c) + ", char = " + str(self.char) + ", field = F" + str(self.q))
            elif self.a == 0 and self.b > 0 and self.c == 0:
                return str("y^2 = x^3 + " + str(self.b) + "x" + "char = " + str(self.char) + ", field = F" + str(self.q))
            else:
                return str("y^2 = x^3" + ", char = " + str(self.char) + ", field = F" + str(self.q))

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

# Helper method to find the inverses of every number mod p. This is used when adding
# two points over a finite field with char p.
def findInverses(p):
    n = 2
    inverses = [0, 1]
    while n < p:
        i = 2
        while i < p:
            if i*n % p == 1:
                inverses.append(i)
            i+=1
        n+=1

    return inverses

# Helper method for factor to tell if an integer is prime.
def isPrime(n):
    if not isinstance(n, int) or isinstance(n, float):
        return "Please pass an int or a float."

    n = int(n)
    l = int(n**(1/2) + 1)
    i = 2
    while i < l:
        if n % i == 0:
            return 0
        i+=1

    return 1


# Helper method to get the prime factors of n. Used in baby step, giant step for point order.
def factor(n):
    if not isinstance(n, int) or isinstance(n, float):
        return "Please pass an int or a float."

    n = int(n)
    if n == 1:
        return [1]

    factors = [1]
    m = int(n**(1/2) + 1)
    if (n % 2) == 0:
        factors.append(2)
    j = 3
    while j <= m:
        if isPrime(j) and (n % j) == 0:
            factors.append(j)
        j+=1

    j = int(n**(1/2) + 1)
    while j < int(n/2 + 1):
        for p in factors:
            if int(p * j) == n:
                factors.append(j)
        j+=1

    if isPrime(n):
        factors.append(n)

    return factors
