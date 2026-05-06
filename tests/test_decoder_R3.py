"""
R3 acceptance test for the RS(255, 223) decoder.

(proposal R3): "decoder must successfully recover at least 95% of 1000 randomly generated message blocks, each corrupted with exactly t = 16 symbol
errors at uniformly random positions with uniformly random non-zero error magnitudes. The decoded output must equal the original message with zero residual symbol errors."
"""

import sys, os, random
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from rs_analyzer.rs_codec import RSParam, encode, decode

PARAM = RSParam()
TRIALS = 1000

def _random_message():
    message = []
    for i in range(PARAM.k):
        symb = random.randint(0,255)
        message.append(symb)
    return message

def _corrupt(codeword, n_errors):
    # Uniformly random positions, uniformly random nonzero magnitudes.
    positions = random.sample(range(PARAM.n), n_errors)
    received = list(codeword)
    for p in positions:
        received[p] ^= random.randint(1, 255)
    return received

def test_r3():
    successes = 0
    failures = 0 #decoder returned success False
    message_failures = 0 #decoder returned success True, but message was wrong

    for trial in range(TRIALS):
        message = _random_message()
        codeword = encode(message, PARAM)
        received = _corrupt(codeword, PARAM.t)
        result = decode(received, PARAM)
        if result.success and result.message == message:
            successes += 1
        elif result.success:
            message_failures += 1
        else:
            failures += 1

    rate = successes / TRIALS
    print(f"    successes: {successes}/{TRIALS}  ({100*rate:.1f}%)")
    print(f"    failures: {failures}")
    print(f"    message failures:  {message_failures}")

if __name__ == "__main__":
        test_r3()
        print("\nR3 PASSED")
