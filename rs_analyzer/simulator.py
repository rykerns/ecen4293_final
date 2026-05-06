import numpy as np
import math, random

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from rs_analyzer.rs_codec import RSParam, encode, decode

# --- Bit and symbol packing ---

"""
Each GF(2^8) symbol is 8 bits, A codeword of n = 255 symbols is therefore 255 * 8 = 2040 bits. 
We need to pass MSB first in packing (though from my research it doesnt seem to matter as long as its consistent between packing and unpacking methods)
"""

def sym_to_bit(symbols):
    """Unpack list of GF symbols into a 1D array. Output lenght is 8* len(symb) MSB first"""
    bits = np.zeros(8 * len(symbols), dtype=np.uint8)
    for i, s in enumerate(symbols):
        for b in range(8):
            bits[8 * i + b] = (s >> (7 - b)) & 1
    return bits

def bit_to_sym(bits):
    """Pack a 1D array into a list of GF symbols. Input length must be a multiple of 8"""

    if len(bits) % 8 != 0:
        raise ValueError(f"bit array length {len(bits)} not a multiple of 8")
    n_symbols = len(bits) // 8
    symbols = [0] * n_symbols
    for i in range(n_symbols):
        s = 0
        for b in range(8):
            s = (s << 1) | int(bits[8 * i + b])
        symbols[i] = s
    return symbols

# --- BSPK + AWGN Channel Simulator ---

"""
BPSK maps {0, 1} -> {+1, -1}. The transmitted signal has unit energy per channel symbol (per bit, since BPSK is one bit per symbol).

Noise variance, accounting for code rate R = k/n: the noise per channel is sigma^2 = N0 / 2 = 1 / (2 * R * (Eb/N0)_linear). The R factor is
critical - omitting it gives a coding gain off by 10*log10(R) ~= 0.6 dB for RS(255, 223)

Hard-decision demodulation: receiver decides 0 if received sample > 0, else 1. This is optimal for uncoded BPSK and is considered standard for an RS decoder.
"""

def bpsk_awgn(bits, ebn0_db, rate, rng):
    """Transmit bits over BPSK + AWGN at given Eb/N0 (dB), return received bits."""
    symbols = 1.0 - 2.0 * bits.astype(np.float64)

    # Noise variance per channel use, with code-rate correction.
    ebn0_linear = 10 ** (ebn0_db / 10.0)
    sigma2 = 1.0 / (2.0 * rate * ebn0_linear)
    sigma = math.sqrt(sigma2)

    noise = rng.normal(0.0, sigma, size=len(symbols))
    received = symbols + noise

    # Hard decision: positive -> 0, negative -> 1
    return (received < 0).astype(np.uint8)

def monte_carlo_coded_ber(params, ebn0_db, max_codewords, max_errors, seed = None):
    """Measure post-decoding bit error rate (BER) for an RS encoded BSPK link"""
    rate = params.k / params.n
    rng = np.random.default_rng(seed) # generative AI suggested adding in a seed for reproducibility

    n_bit_errors = 0
    n_bits = 0
    n_decode_failures = 0
    n_codewords = 0

    for l in range(max_codewords):
        # Generate random message, encode, modulate, transmit, decode
        message = []
        for i in range(params.k):
            symb = rng.integers(0,256) # turns out rng in numpy is upper bound exclusive, unlike in the random library
            message.append(symb)
        codeword = encode(message, params)

        tx_bits = sym_to_bit(codeword)
        rx_bits = bpsk_awgn(tx_bits, ebn0_db, rate, rng)
        received = bit_to_sym(rx_bits)

        result = decode(received, params)
        if not result.success:
            n_decode_failures += 1
            # On decode failure, the message field of 'received' is the decoder's best (uncorrected) guess at the message symbols.
            # Count bit errors against that to keep the BER measurement throwing the codeword away would underreport errors.
            decoded_message = received[params.parity_len:]
        else:
            decoded_message = result.message

        # Compare decoded message bits to original message bits.
        original_bits = sym_to_bit(message)
        decoded_bits = sym_to_bit(decoded_message)
        n_bit_errors += int(np.sum(original_bits != decoded_bits))
        n_bits += len(original_bits)
        n_codewords += 1

        if n_bit_errors >= max_errors:
            break

    return {
        "ebn0_db": ebn0_db,
        "n_bit_errors": n_bit_errors,
        "n_bits": n_bits,
        "n_codewords": n_codewords,
        "n_decode_failures": n_decode_failures,
        "ber": n_bit_errors / n_bits if n_bits > 0 else 0.0,
    }

def monte_carlo_uncoded_ber(ebn0_db, n_bits, seed=None):
    """Measure raw BPSK + AWGN bit-error rate (no coding) for comparison.

    Should agree with Q(sqrt(2 * Eb/N0)) within Monte Carlo noise. Used as a sanity check on the channel implementation to get data to test against

    """
    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, size=n_bits, dtype=np.uint8)
    rx_bits = bpsk_awgn(bits, ebn0_db, rate=1.0, rng=rng)
    n_bit_errors = int(np.sum(bits != rx_bits))
    return {
        "ebn0_db": ebn0_db,
        "n_bit_errors": n_bit_errors,
        "n_bits": n_bits,
        "ber": n_bit_errors / n_bits,
    }