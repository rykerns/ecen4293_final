"""
Acts as the RS encoder, decoder is later (when we implement berlekamp-massey/other numerical methods) but needs to share some kind of RS parameter object (class probably) with encoder

Same field conventions as gf256.py

Generator polynomial:
    g(x) = prod_{i=1}^{2t} (x - alpha^i)
(the first consecutive root is alpha^1 at b = 1 convention).

messages and codewords are lists of ints in [0,255]; a message has length k; a codeword has length n=k+2t

Codeword layout is parity in the LOW positions:
    codeword[0 : 2t] = parity symbols
    codeword[2t : n] = original message symbols

Polynomial-wise, with poly[i] = coefficient of x^i:
    c(x) = x^{2t} * m(x) + r(x)
where r(x) = x^{2t} * m(x) mod g(x).
"""

from rs_analyzer import gf256
from dataclasses import dataclass, field
from typing import Optional

@dataclass
class RSParam:
    """
    Parameters for an RS(n, k) code over GF(2^8).
    Define: n = 255 symbols and k = 223 message symbols, capable of correcting up to t=16 symbol errors per codeword, k=n-2t
    Needs to pass:
     - error correction capacity: t = (n-k)/2 symbols
     - a constructed list of g(x)=prod_{i=1}^{2t} (x - alpha^i) over GF(2^8)
    """

    # Parameters for the RSParam class for RS(n,k) code
    n: int = 255
    k: int = 223
    generator: list = field(init = False, repr = False)

    # we must construct this such that k<0 and n<k (this is true for this RS encoding, but now we can vary the implementation and try different numbers)
    def __post_init__(self):
        if not (0 < self.k < self.n <= 255):
            raise ValueError(f"require 0 < k < n <= 255 -- recieved n = {self.n}, k = {self.k}")
        # need n-k even for clean 2t pairity symbols
        if (self.n - self.k) % 2 != 0:
            raise ValueError(f"n-k must be even -- recieved n-k = {self.n - self.k}")
        self.generator = _build_generator(self.t)

    @property
    def t(self) -> int:
        # Total error correction capacity computed from n and k
        cap = (self.n - self.k) // 2
        return cap
    
    @property
    def parity_len(self) -> int:
        # No. of pairity symbols: 2t = n - k
        res = self.n - self.k
        return res
    
def _build_generator(t: int) -> list:
    """
    Construct g(x) = prod_{i=1}^{2t} (x - a^i) over GF(2^8).

    In GF(2^8), (x - a) = (x + a), so each factor is [a, 1] in little-endian polynomial convention.
    """

    g = [1]  # start with g(x) = 1
    for i in range(1, 2 * t + 1):
        alpha_i = gf256.power(2, i)  # alpha = 2
        g = gf256.poly_mul(g, [alpha_i, 1])
    return g

# --- Encoder part ---

def encode(message, params: RSParam) -> list:
    """
    Encodes a length-k message into a length-n codeword

    # Polynomial view: m(x) with message[0] as the constant term
    # x^{2t} * m(x) places the message into the high-order coeff, leaving the low 2t positions zero for the parity to occupy.
    """
    _validate_message(message, params)

    # place message into the higher order positions and append the pairity bits
    shifted = [0] * params.parity_len + list(message)

    # remainder r(x) = shifted(x) mod g(x)
    _, r = gf256.poly_divmod(shifted, params.generator)

    #build codeword: start from shifted, XOR in the remainder at the low 2t positions
    codeword = list(shifted)
    for i, coeff in enumerate(r):
        codeword[i] ^= coeff
    return codeword

def _validate_message(message, params: RSParam) -> None:
    """Throw errors if either message length does not match k OR if the message is outside the GF(2^8) range"""

    if len(message) != params.k:
        raise ValueError(f"Message length {len(message)} does not match expected k: {params.k} ")
    for i, seq in enumerate(message): 
        if not (0 <= seq <= 255):
            raise ValueError(f"message[{i}] = {seq} is outside the GF(2^8) range [0,255]")

def syndromes(recieved, params: RSParam) -> list:
    """
    Compute the 2t syndromes S_i = received(alpha^i) for i = 1..2t.
 
    A valid codeword yields all-zero syndromes because it is divisibleby g(x), whose roots are exactly alpha^1, ... alpha^{2t}.
    """
    s_list = []
    for i in range(1, params.parity_len + 1):
        alpha_i = gf256.power(2, i)
        s_i = gf256.poly_eval(recieved, alpha_i)
        s_list.append(s_i)
    return s_list

def is_codeword(word, params:RSParam):
    """True iff `word` has length n and all 2t syndromes are zero."""
    if len(word) != params.n:
        return False
    for s in syndromes(word, params):
        if s != 0:
            return False
    return True

# --- Solver Functions ---

def berlekamp_massey(syndromes_list, params:RSParam) -> list:
    """
    Recover the error locator polynomial Lambda(x) frp, the syndrome squence 

    Given syndromes S_1, S_2, ..., S_{2t}, returns the minimum-degree polynomial Lambda(x) such that the linear recurrence:
        sum_{j=0}^{L} Lambda[j] * S[i-j] = 0
    holds for i > L, where L = deg(Lambda)

    In little-endian polynomial convention, Lambda[0] is the constant term and equals 1 (the polynomial is monic when reversed)
    The roots of Lambda are at alpha^{-i_j} for each error position i_j

    If the input syndromes come from a received word with at most t errors, deg(Lambda) equals the number of errors
    If more than t errors occurred, Lambda may have degree > t or may fail to factor
    """

    if len(syndromes_list) != params.parity_len:
        raise ValueError(f"expected {params.parity_len} syndromes, got {len(syndromes_list)}")
    
    S = syndromes_list
    n_syndrome = len(S)

    # If all syndromes are zero, the received word is already a codeword and the trivial locator Lambda(x) = 1 is correct.
    if all(s == 0 for s in S):
        return [1]
    
    Lambda = [1]# current locator polynomial (Lambda(0) = 1)
    B = [1] # previous "best" polynomial before last length update
    L = 0  # current LFSR length
    m = 1  # iterations since last length update
    b = 1  # discrepancy at last length update

    # iterate over syndromes, 0 indexed so S[i] is the current syndrome being matched
    for i in range(n_syndrome):
        # check for discrepancy delta = S[i] + sum_{j=1...L} Lambda[j] * S[i - j].
        delta = S[i]
        for j in range(1, L+1):
            if j< len(Lambda) and Lambda[j] != 0:
                delta ^= gf256.mul(Lambda[j], S[i-j])
        
        # if delta is 0, then Lambda already predicts S, so we just advance m
        if delta == 0:
            m += 1
            continue
        
        # Correction term is (delta / b) * x^m * B, XOR into Lambda
        coef = gf256.div(delta, b)
        # Build correction = coef * x^m * B, then XOR into Lambda
        correction_len = m + len(B)

        if correction_len > len(Lambda):
            Lambda = Lambda + [0] * (correction_len - len(Lambda))
        new_Lambda = list(Lambda)
        for j in range(len(B)):
            if B[j] != 0:
                new_Lambda[j + m] ^= gf256.mul(coef, B[j])

        # Length update condition: 2*L <= i
        if 2 * L <= i:
            T = Lambda # save current Lambda before overwrite
            Lambda = new_Lambda
            B = T # B becomes pre-update Lambda
            L = i + 1 - L # new length
            b = delta
            m = 1
        else:
            Lambda = new_Lambda
            m += 1

    # Strip trailing zeros from the result
    while len(Lambda) > 1 and Lambda[-1] == 0:
        Lambda.pop()
    return Lambda

def chien_search(Lambda, params:RSParam):
    """
    Find error positions by locating roots of Lambda over GF(2^8)

    For each position p in [0, n), test whether alpha^(-p) is a root of Lambda. Each root corresponds to an error at position p
    """
    deg = len(Lambda) - 1
    if deg == 0:
        return [] # Lambda(x) = 1: no errors to locate

    positions = []
    for p in range(params.n):
        x_inv = gf256.power(2, -p)  # alpha^(-p)
        if gf256.poly_eval(Lambda, x_inv) == 0:
            positions.append(p)

    if len(positions) != deg:
        return None
    return positions


def forney(Lambda, error_pos, syndromes_list, params: RSParam):
    """
    Compute error magnitudes at the given positions via Forney's formula.

    With b = 1 convention:
        e_j = Omega(alpha^{-i_j}) / Lambda'(alpha^{-i_j})
    where Omega(x) = [S(x) * Lambda(x)] mod x^{2t}.

    syndromes_list is in little-endian order (index j holds S_{j+1}), so it serves as the coefficient list of S(x).
    """

    # build error evaluator polynomial
    Omega_full = gf256.poly_mul(syndromes_list, Lambda) #unreduced product: If S has degree 2t-1 and Lambda has degree ν <= t, the product has degree up to (2t-1) + ν, which is bigger than 2t−1. So this has more coeff than we want
    Omega = Omega_full[:params.parity_len] #in little endian form the coeff x^i is at i, so we keep the indices [0,2t)

    Lambda_prime = gf256.poly_deriv(Lambda)

    magnitudes = []
    for p in error_pos:
        x_inv = gf256.power(2, -p) # alpha^{-i_j}
        num = gf256.poly_eval(Omega, x_inv)
        den = gf256.poly_eval(Lambda_prime, x_inv)
        if den == 0:
            # Lambda'(x) vanishing at a root of Lambda means Lambda has a repeated root, which cannot arise from a valid error pattern.
            # Treat as a decode failure rather than divide by zero.
            raise ValueError(f"Forney: Lambda' vanishes at alpha^(-{p}); locator does not factor cleanly")
        magnitudes.append(gf256.div(num, den))
    return magnitudes

# --- Decoder ---

@dataclass
class DecodeResult:
    """
    Outcome of a decoder attempt
    
    success=True with n_errors_corrected=k means the decoder applied k corrections and the result is internally consistent (zero residual syndromes). 
    
    It does NOT prove the result equals the original message: a >t-error pattern can occasionally produce a Lambda that factors and
    a correction that zeros syndromes but lands on the wrong codeword. That case is undetectable from inside the decoder
    """
    message: Optional[list]
    success: bool
    n_errors_corrected: int = 0
    failure_reason: Optional[str] = None

def decode(received, params:RSParam):
    """
    Full RS syndrome decoder Runs: syndromes -> Berlekamp-Massey -> Chien -> Forney -> XOR corrections
    and returns a DecodeResult distinguishing the failure modes from a clean decode.

    Go top down from fastest failure to slowest (likely saves compute)
    """

    if len(received) != params.n:
        raise ValueError(f"received length {len(received)} does not match n = {params.n}")
    
    # zero-syndrome is already a code word
    S = syndromes(received, params)
    all_zero = True
    for s in S:
        if s != 0:
            all_zero = False
            break
    if all_zero:
        message = list(received[params.parity_len:])
        return DecodeResult(message=message, success = True, n_errors_corrected=0)
    
    # if not, run the error locator
    Lambda = berlekamp_massey(S, params)

    # If BM returns degree > t, more than t errors occurred. No point running Chien - the locator can't correspond to a correctable pattern.
    if len(Lambda) -1 > params.t:
        return DecodeResult(message = None, success = False, failure_reason = f"deg(Lambda) = {len(Lambda)-1} exceeds t = {params.t}")
    
    """
    If not, run error positions check with Chien. Chien evaluates Lambda at all 255 nonzero field elements and counts roots. Three things can happen:
     - root count = deg(L): clean factorization, every root corresponds to a real error position. "good ending"
     - root count < deg(L): Lambda doesn't fully factor over GF(2^8). It has roots at infinity or roots in some extension field. The polynomial BM returned isn't the one a real <= t-error pattern would produce
     - root count > deg(L): impossible by the fundamental theorem of algebra (a polynomial has at most deg roots)
    """
    pos = chien_search(Lambda, params)
    if pos is None:
        return DecodeResult(message = None, success = False, failure_reason = "Chien roots != deg(Lambda), locator does not factor >t errors")
    
    # if not then we check forney, unique check for when Lambda has a repeated root so two factors collapsed into one because i_a = i_b. Likely any repeated roots wont get past Chien, but this produces a clean Decode result, rather than a 'divide by zero error'

    try:
        mag = forney(Lambda, pos, S, params)
    except ValueError as ex:
        return DecodeResult(message = None, success = False, failure_reason = f"Forney caught: {e}")
    
    # Then we apply corrections:
    corrected = list(received)
    for p, e in zip(pos, mag):
        corrected[p] ^= e

    # After correction, syndromes must still be zero. If they arent, the pattern was a >t-error case that produced a plausible Lambda but a wrong correction.
    corr_check = syndromes(corrected, params)
    any_nonzero = False
    for s in corr_check:
        if s != 0:
            any_nonzero = True
            break
    if any_nonzero:
        return DecodeResult(message = None, success = False, failure_reason = "Post-correction syndromes contained zero. Edge case failure inside BM or Forney.")
    
    message = corrected[params.parity_len:]
    return DecodeResult(message = message, success = True, n_errors_corrected = len(pos))

