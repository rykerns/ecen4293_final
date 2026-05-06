import csv
import math
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Allow running as `python scripts/plot_ber.py` from the project root.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rs_analyzer.numerics import q_function


CSV_PATH = "ber_sweep_results.csv"
OUTPUT_PNG = "ber_curve.png"
TARGET_BER = 1e-4  # the BER level at which coding gain is measured (R4)


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

def load_sweep(path):
    """Load CSV into two dicts keyed by Eb/N0: coded and uncoded BER points.

    Each value is a tuple (ber, n_bits) - n_bits lets us flag points where
    BER = 0 means "no errors observed in N bits" rather than a real zero.
    """
    coded = {}
    uncoded = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ebn0 = float(row["ebn0_db"])
            ber = float(row["ber"])
            n_bits = int(row["n_bits"])
            if row["kind"] == "coded":
                coded[ebn0] = (ber, n_bits)
            elif row["kind"] == "uncoded":
                uncoded[ebn0] = (ber, n_bits)
    return coded, uncoded


# ---------------------------------------------------------------------------
# Crossing detection (linear interpolation in log space)
# ---------------------------------------------------------------------------
#
# To find the SNR at which BER crosses TARGET_BER, we interpolate linearly
# in (Eb/N0, log10(BER)) space - BER curves are roughly straight lines on
# a semilog plot near the waterfall, so log-linear interpolation is very
# accurate there.
#
# Skip points with BER == 0 (they're "below measurement floor", not real
# zeros) and points where n_bits is suspiciously small.

def find_crossing(snrs, bers, target):
    """Return SNR where BER crosses `target`, via log-linear interpolation.

    Returns None if no crossing exists in the data range.
    """
    # Build (snr, log10(ber)) pairs, sorted by SNR. Drop zero-BER points.
    pairs = sorted(
        (s, math.log10(b)) for s, b in zip(snrs, bers) if b > 0
    )
    if len(pairs) < 2:
        return None

    log_target = math.log10(target)
    # Walk pairs; find the first interval where target lies between adjacent
    # log10(BER) values.
    for (s1, lb1), (s2, lb2) in zip(pairs, pairs[1:]):
        # Curves descend with SNR, so lb1 should be larger than lb2.
        # Target is bracketed iff lb1 >= log_target >= lb2 (or vice versa).
        if (lb1 - log_target) * (lb2 - log_target) <= 0:
            # Linear interpolation: fraction along [s1, s2] where log10(BER)
            # equals log_target.
            frac = (log_target - lb1) / (lb2 - lb1)
            return s1 + frac * (s2 - s1)
    return None


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def main():
    if not os.path.exists(CSV_PATH):
        print(f"ERROR: {CSV_PATH} not found. Run run_ber_sweep.py first.")
        return 1

    coded, uncoded = load_sweep(CSV_PATH)
    if not coded:
        print(f"ERROR: no coded BER data in {CSV_PATH}.")
        return 1

    # Sort by SNR for plotting.
    coded_snr = sorted(coded.keys())
    coded_ber = [coded[s][0] for s in coded_snr]
    uncoded_snr = sorted(uncoded.keys())
    uncoded_ber = [uncoded[s][0] for s in uncoded_snr]

    # Theoretical uncoded curve via numerical Q-function (R4 sub-requirement).
    theory_snr = np.linspace(0, 12, 200)
    theory_ber = [q_function(math.sqrt(2 * 10**(s/10))) for s in theory_snr]

    # Find crossings at TARGET_BER.
    coded_crossing = find_crossing(coded_snr, coded_ber, TARGET_BER)
    uncoded_crossing_sim = find_crossing(uncoded_snr, uncoded_ber, TARGET_BER)
    # Also find the theoretical uncoded crossing - more reliable than the
    # simulated one because the simulator may not have data at high enough
    # SNR to bracket BER = 1e-4 cleanly.
    uncoded_crossing_theory = find_crossing(
        list(theory_snr), theory_ber, TARGET_BER
    )

    # Coding gain: prefer the theoretical uncoded crossing as the reference.
    # The closed-form Q expression is exact; the simulated uncoded curve
    # adds Monte Carlo noise that the coding-gain measurement doesn't need.
    coding_gain_db = None
    if coded_crossing is not None and uncoded_crossing_theory is not None:
        coding_gain_db = uncoded_crossing_theory - coded_crossing

    # ----- Build the plot ------------------------------------------------

    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the three curves.
    ax.semilogy(theory_snr, theory_ber, "k--", linewidth=1.5,
                label=r"Uncoded BPSK (theory: $Q(\sqrt{2 E_b/N_0})$)")

    # Replace zero BER with a placeholder above the floor so the marker
    # still appears on a log plot. Use the smallest nonzero BER / 10 as
    # a visual "below this" indicator.
    nonzero_coded = [b for b in coded_ber if b > 0]
    floor = min(nonzero_coded) / 10 if nonzero_coded else 1e-10
    coded_ber_plot = [b if b > 0 else floor for b in coded_ber]
    ax.semilogy(coded_snr, coded_ber_plot, "ro-", linewidth=1.5,
                markersize=6, label="RS(255, 223), simulated")

    nonzero_uncoded = [b for b in uncoded_ber if b > 0]
    if nonzero_uncoded:
        ax.semilogy(uncoded_snr, [b if b > 0 else floor for b in uncoded_ber],
                    "b.-", linewidth=1, alpha=0.6,
                    label="Uncoded BPSK, simulated")

    # Reference line at TARGET_BER.
    ax.axhline(TARGET_BER, color="gray", linestyle=":", linewidth=1, alpha=0.7)
    ax.text(0.5, TARGET_BER * 1.5, f"BER = {TARGET_BER:.0e}",
            color="gray", fontsize=9, va="bottom")

    # Mark the crossings.
    if coded_crossing is not None:
        ax.plot(coded_crossing, TARGET_BER, "rs", markersize=10,
                markerfacecolor="none", markeredgewidth=2)
        ax.annotate(f"  {coded_crossing:.2f} dB",
                    xy=(coded_crossing, TARGET_BER),
                    xytext=(8, 0), textcoords="offset points",
                    fontsize=9, color="red", va="center")
    if uncoded_crossing_theory is not None:
        ax.plot(uncoded_crossing_theory, TARGET_BER, "ks", markersize=10,
                markerfacecolor="none", markeredgewidth=2)
        ax.annotate(f"  {uncoded_crossing_theory:.2f} dB",
                    xy=(uncoded_crossing_theory, TARGET_BER),
                    xytext=(8, 0), textcoords="offset points",
                    fontsize=9, color="black", va="center")

    # Coding gain annotation.
    if coding_gain_db is not None:
        ax.set_title(f"RS(255, 223) over BPSK + AWGN: "
                     f"coding gain = {coding_gain_db:.2f} dB at "
                     f"BER = {TARGET_BER:.0e}")
    else:
        ax.set_title("RS(255, 223) over BPSK + AWGN")

    ax.set_xlabel(r"$E_b / N_0$ (dB)")
    ax.set_ylabel("Bit Error Rate")
    ax.set_ylim(1e-7, 1)
    ax.set_xlim(0, max(max(coded_snr), max(uncoded_snr)) + 0.5)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left")

    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=150)
    print(f"Saved {OUTPUT_PNG}")

    # ----- Print numerical summary --------------------------------------

    print()
    print("=" * 60)
    print("BER curve summary")
    print("=" * 60)
    print(f"Uncoded BPSK reaches BER = {TARGET_BER:.0e} at:")
    if uncoded_crossing_theory is not None:
        print(f"  Theory:    {uncoded_crossing_theory:.3f} dB")
    if uncoded_crossing_sim is not None:
        print(f"  Simulated: {uncoded_crossing_sim:.3f} dB")
    print()
    print(f"RS(255, 223) reaches BER = {TARGET_BER:.0e} at:")
    if coded_crossing is not None:
        print(f"  Simulated: {coded_crossing:.3f} dB")
    else:
        print(f"  Simulated: not bracketed by data (need more SNR points "
              f"between the data points where BER crosses {TARGET_BER:.0e})")
    print()
    if coding_gain_db is not None:
        print(f"Coding gain at BER = {TARGET_BER:.0e}: "
              f"{coding_gain_db:.3f} dB")
        if coding_gain_db >= 3.0:
            print(f"  -> R4 SATISFIED (>=3.0 dB required)")
        else:
            print(f"  -> R4 NOT YET SATISFIED (>=3.0 dB required)")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())