import csv
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rs_analyzer.rs_codec import RSParam
from rs_analyzer.simulator import monte_carlo_coded_ber, monte_carlo_uncoded_ber

SNR_POINTS_DB = [
    # Below the waterfall - decoder fails almost always, sparse sampling.
    1.0, 2.0, 3.0, 4.0, 4.5,
    # Waterfall region - dense sampling, this is where the curve is steep.
    5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5,
    # Above the waterfall - decoder succeeds almost always, BER very small.
    7.0, 7.5, 8.0,
]

MAX_CODEWORDS = 50000
MAX_ERRORS = 200  # stop accumulating once 200 bit errors observed

OUTPUT_CSV = "ber_sweep_results.csv"
PARAMS = RSParam()

def main():
    print(f"Sweep configuration:")
    print(f"    SNR points: {SNR_POINTS_DB}")
    print(f"    max_codewords: {MAX_CODEWORDS}")
    print(f"    max_errors: {MAX_ERRORS}")
    print(f"    output CSV: {OUTPUT_CSV}\n")

    # Open CSV in append mode if it already exists; write header only if the file is new. Provided by generative AI
    file_exists = os.path.exists(OUTPUT_CSV)
    with open(OUTPUT_CSV, "a", newline="") as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow([
                "kind", "ebn0_db", "n_bit_errors", "n_bits",
                "n_codewords", "n_decode_failures", "ber", "elapsed_sec",
            ])
            f.flush()

        for i, ebn0 in enumerate(SNR_POINTS_DB):
            t_start = time.time()

            # Coded run.
            coded = monte_carlo_coded_ber(
                PARAMS, ebn0_db=ebn0,
                max_codewords=MAX_CODEWORDS,
                max_errors=MAX_ERRORS,
                seed=42 + i,  # reproducible, but different per SNR point
            )
            elapsed = time.time() - t_start

            writer.writerow([
                "coded", ebn0, coded["n_bit_errors"], coded["n_bits"],
                coded["n_codewords"], coded["n_decode_failures"],
                coded["ber"], f"{elapsed:.1f}",
            ])
            f.flush()  # ensures the row hits disk now, not at script exit

            print(f"[{i+1}/{len(SNR_POINTS_DB)}] coded @ {ebn0:.2f} dB: "
                  f"BER = {coded['ber']:.3e}  "
                  f"({coded['n_codewords']} cw, "
                  f"{coded['n_decode_failures']} fails, {elapsed:.1f}s)")

            # Uncoded reference. Cheap - just BPSK + AWGN, no codec.
            t_start = time.time()
            uncoded = monte_carlo_uncoded_ber(
                ebn0, n_bits=2_000_000, seed=42 + i,
            )
            elapsed = time.time() - t_start

            writer.writerow([
                "uncoded", ebn0, uncoded["n_bit_errors"], uncoded["n_bits"],
                "", "",  uncoded["ber"], f"{elapsed:.1f}"])
            f.flush()

            print(f"           uncoded @ {ebn0:.2f} dB: "f"BER = {uncoded['ber']:.3e}  ({elapsed:.1f}s)")
    print()
    print(f"Sweep complete. Results in {OUTPUT_CSV}.")

if __name__ == "__main__":
    main()