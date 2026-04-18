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

# --- Scalar operations ---

def add(a, b):
    """GF(2^8) addition"""
    return a ^ b

def sub(a, b):
    """GF(2^8) subtraction; Identical to addition in characteristic 2. Decided to keep both in since the add/sub inverse rule is only true for certian generator alpha."""
    return a ^ b

def mul(a, b):
    """GF(2^8) multiplication via log/antilog tables."""
    if a == 0 or b == 0:
        return 0
    return EXP[LOG[a] + LOG[b]]

def div(a, b):
    """GF(2^8) division; throw error if b == 0."""
    if b == 0:
        raise ZeroDivisionError("division by zero in GF(2^8)")
    if a == 0:
        return 0
    return EXP[LOG[a] - LOG[b] + ORDER]

def inv(a):
    """GF(2^8) multiplicative inverse; throw error if a == 0."""
    if a == 0:
        raise ZeroDivisionError("zero has no multiplicative inverse in GF(2^8)")
    return EXP[ORDER - LOG[a]]

def power(a, n):
    """
    GF(2^8) exponentiation a^n for any integer n.

    0^0 = 1. 0^n for n > 0 is 0. 0^n for n < 0 raises.
    """
    if a == 0:
        if n == 0:
            return 1
        if n < 0:
            raise ZeroDivisionError("zero raised to a negative power")
        return 0
    return EXP[(LOG[a] * n) % ORDER]

# --- Polynomial Operations ---

def _strip(p):
    """Strip trailing zero coefficients, keeping at least [0]."""
    out = list(p)
    while len(out) > 1 and out[-1] == 0:
        out.pop() # trim high-degree zeros
    return out

def poly_add(p, q):
    """Add two polynomials over GF(2^8)."""
    n = max(len(p), len(q)) #length of result is the longer operand
    result = [0] * n
    for i in range(len(p)):
        result[i] ^= p[i] #XOR  in p coeff
    for i in range(len(q)):
        result[i] ^= q[i] #XOR in q coeff
    return _strip(result)

def poly_scale(p, c):
    """Multiply polynomial p by scalar c. in a field there are no zero divisors, so if the input was already stripped and c =/= 0, the output is automatically stripped too."""
    if c == 0:
        return [0] #anything * 0 is the zero poly
    if c == 1:
        return list(p) #anything * 1 is a copy
    lc = LOG[c] # precompute LOG[c] once, we could call mul
    result = []
    for coef in p: # loop through coeffs of poly
        if coef == 0: #check for non-zero elements
            result.append(0)
        else:
            result.append(EXP[LOG[coef] + lc]) #LOG[coeff] is the discrete log of the coefficient: which power of alpha = coeff
    return result

def poly_mul(p, q):
    """Multiply two polynomials over GF(2^8); is essentially convolution. Multiplying a degree m polynomial by a degree n polynomial gives degree m+n, so the result has m+n+1 coefficients"""
    if not p or not q: 
        return [0]
    result = [0] * (len(p) + len(q) - 1) #degree of resulting polynomial
    for i, a in enumerate(p):
        if a == 0:
            continue #skip the zero terms in p
        la = LOG[a]
        for j, b in enumerate(q):
            if b == 0:
                continue #skip the zero terms in q
            result[i + j] ^= EXP[la + LOG[b]] #accumulate a*b into index i+j, need this since multiple (i,j) pairs land at the same k, so their contributions XOR into that slot
    return _strip(result)

def poly_eval(p, x):
    """Evaluate polynomial p at element x using Horner's method: instead of computing powers of x^n seperatly, we nest each polynomial, then repeatedly "multiply by x and add the next coefficient."""
    if not p:
        return 0 # skip zeros
    result = p[-1] #this is the highest degree coeff since  p = [1, 0 ,2, 3] = 1+2x^2+3x^3 so p[-1] = 3
    for i in range(len(p) - 2, -1, -1): #go DOWN towards i=0
        result = mul(result, x) ^ p[i] # stack the results, "^" acts as addition in the field
    return result

def poly_divmod(num, den):
    """
    Polynomial long division over GF(2^8).

    Returns (quotient, remainder) such that num == poly_add(poly_mul(quotient, den), remainder) with deg(remainder) < deg(den).
    """
    den = _strip(den) #ensure den[-1] is the leading coeff of the denominator
    if den == [0]: #throw error if zero even after strip
        raise ZeroDivisionError("division by zero polynomial") 
    num = list(num) 
    if not num:
        return [0], [0] #easy: zero case
    deg_den = len(den) - 1
    if len(num) - 1 < deg_den:
        return [0], _strip(num) # easy: if the num is already smaller than den: q=0 and r=num 

    inv_lead = inv(den[-1]) #compute leading coeff of  1/den(max coeff)
    quotient = [0] * (len(num) - deg_den) # space for deg(num)-deg(den)+1 coefficients
    while len(num) - 1 >= deg_den:
        leading = num[-1]
        if leading == 0:
            num.pop() #skip zeros, and remove from num
            continue
        q_deg = len(num) - 1 - deg_den #check what power of x this quotient term has
        q_coef = mul(leading, inv_lead) # 1/d
        quotient[q_deg] = q_coef
        for i in range(len(den)):
            num[q_deg + i] ^= mul(q_coef, den[i]) # num -= q_coef * x^q_deg * den
        num.pop()  # top coefficient is now zero by construction
    return _strip(quotient), _strip(num) if num else [0]

def poly_deriv(p):
    """
    Formal derivative of polynomial p over GF(2^8).
 
    In characteristic 2, d/dx(x^i) = x^(i-1) if i is odd, else 0.

    In GF(p) for odd p you'd multiply each coefficient by i mod p; in GF(2ⁿ) the factor is always 0 or 1
    """
    if len(p) <= 1: 
        return [0] # derivative of constant is zero
    result = [0] * (len(p) - 1) # one fewer coeff
    for i in range(1, len(p)): #skip constant in all other cases
        if i & 1: # is i odd?
            result[i - 1] = p[i] #if so keep and shift by one position
    return _strip(result)

