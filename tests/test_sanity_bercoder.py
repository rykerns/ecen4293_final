import numpy as np
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rs_analyzer.simulator import (bit_to_sym, sym_to_bit, monte_carlo_uncoded_ber, monte_carlo_coded_ber,)
from rs_analyzer.rs_codec import RSParam
from rs_analyzer.numerics import q_function
import math

# Check 1: bit packing roundtrip. If this fails, every BER is junk
test = [0, 1, 42, 200, 255]
assert bit_to_sym(sym_to_bit(test)) == test
print("bit packing roundtrip OK")

# Check 2: uncoded BPSK matches Q(sqrt(2 Eb/N0)) within statistical noise. At Eb/N0 = 5 dB, theory says BER ~= 5.95e-3. With 1e6 bits we expect
# about 6000 errors, so MC noise is roughly 1/sqrt(6000) ~= 1.3% relative.
ebn0 = 5.0
result = monte_carlo_uncoded_ber(ebn0, n_bits=1_000_000, seed=0)
theory = q_function(math.sqrt(2 * 10**(ebn0/10)))
print(f"uncoded @ {ebn0} dB: simulated = {result['ber']:.4e}, theory = {theory:.4e}")
# These should agree to ~3 significant figures.

# Check 3: coded BER at one SNR point, small budget. At 4 dB the coded
# curve should steeply descending; ideally decode_failures > 0 but successes too. This shakes out any decoder/simulator mismatch.
for ebn0 in [4.0, 5.0, 5.5, 6.0, 6.5]:
    result = monte_carlo_coded_ber(RSParam(), ebn0_db=ebn0, max_codewords=500, max_errors=5000, seed=0)
    pct_fail = 100 * result['n_decode_failures'] / result['n_codewords']
    print(f"@ {ebn0:.1f} dB: BER = {result['ber']:.3e}, "f"decode failures = {pct_fail:.0f}%")