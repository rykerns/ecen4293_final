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

class RSParam:
    """
    Parameters for an RS(n, k) code over GF(2^8).
    Define: n = 255 symbols and k = 223 message symbols, capable of correcting up to t=16 symbol errors per codeword, k=n-2t
    Needs to pass:
     - error correction capacity: t = (n-k)/2 symbols
     - a constructed list of g(x)=prod_{i=1}^{2t} (x - alpha^i) over GF(2^8); In GF(2^8), (x - a) = (x + a), so each factor is always in [a, 1]
    """


# --- Encoder part ---

def encode(message, params: RSParam) -> list:
    """
    Encodes a length-k message into a length-n codeword

    # Polynomial view: m(x) with message[0] as the constant term
    # x^{2t} * m(x) places the message into the high-order coeff, leaving the low 2t positions zero for the parity to occupy.
    """

