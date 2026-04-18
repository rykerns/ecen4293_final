"""
GF(2^8) field arithmetic and polynomial operations

Implements the finite field GF(2^8) using the primitive polynomial p(x) = x^8 + x^4 + x^3 + x^2 + 1 (= 0x11D)
with the primitive element alpha = 0x02, which is a generator of the set.

Every non-zero element of the field can be written as alpha^i for some unique i in [0,254]

-Field elements will be python ints in [0,255]
-Polynomials are lists of ints with index 0 = constant term (in little endian): [3, 2, 1]  denotes  3 + 2*x + 1*x^2
-The zero polynomial is [0], functions should strip trailing zero coeff.

"""

PRIM_POLY = 0x11D
ORDER = 255

# EXP[i] = alpha^i. Length is 2*ORDER so that mul can skip the modulo step and any valid sum log[a] + log[b] lies in [0,2*254]
EXP = [0] * (2 * ORDER)

# LOG[a] = discrete log of a in base alpha, for a in [1, 255].
# LOG[0] is undefined
LOG = [0] * 256

def _build_tables():
    """Populate EXP and LOG iterating through alpha=2 through all of GF(2^8)"""
    x = 1
    for i in range(ORDER):
        EXP[i] = x
        LOG[x] = i
        # Check overflow and wrap
        x <<= 1 #shift x
        if x & 0x100: #if x has grown past 255 reduce modulo PRIM POLY
            x ^= PRIM_POLY
    # Wrap EXP so EXP[i] for i in [255, 2*255-1] mirrors EXP[0..254]
    for i in range(ORDER, 2 *ORDER):
        EXP[i] = EXP[i - ORDER]

_build_tables()


"""
Then we have a bunch of field operations for the polynomial arithmatic

Scalar helpers: add, sub, mul, div, inv, pow
 - add and sub are essentially XOR since -1 = 1, so addition is its own inverse
 - pow: exponents live in mod255 cyclic group
Polynomial helpers: poly_add, poly_mul, poly_div, poly_dx

"""
