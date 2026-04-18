"""
Test suite for gf256 module for R1

!!DISCLAIMER!! Suite was made by generative AI, I provided it with a list of about 20 tests I wanted to run and the gf256.py module, 
    it suggested some more tests to verify the module, and the total came out to 58 tests
"""
import sys, os, traceback
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rs_analyzer import gf256


# ---------------------------------------------------------------------------
# Log / antilog table structure
# ---------------------------------------------------------------------------
 
def test_alpha_zero_is_one():
    assert gf256.EXP[0] == 1
 
 
def test_alpha_255_is_one():
    # alpha^255 = 1 since |GF(2^8)^*| = 255
    assert gf256.EXP[255] == 1
 
 
def test_log_exp_roundtrip():
    for i in range(255):
        assert gf256.LOG[gf256.EXP[i]] == i
 
 
def test_exp_log_roundtrip():
    for a in range(1, 256):
        assert gf256.EXP[gf256.LOG[a]] == a
 
 
def test_exp_covers_nonzero_field():
    # Every nonzero element must appear exactly once in EXP[0:255].
    assert set(gf256.EXP[:255]) == set(range(1, 256))
 
 
def test_alpha_is_two():
    assert gf256.EXP[1] == 2
 
 
def test_alpha_to_8_matches_primitive_poly():
    # Under p(x) = x^8 + x^4 + x^3 + x^2 + 1, reducing x^8 gives
    # x^4 + x^3 + x^2 + 1 = 0x1D.
    assert gf256.EXP[8] == 0x1D
 
 
# ---------------------------------------------------------------------------
# Addition / subtraction
# ---------------------------------------------------------------------------
 
def test_additive_identity():
    for a in range(256):
        assert gf256.add(a, 0) == a
        assert gf256.add(0, a) == a
 
 
def test_add_self_inverse():
    for a in range(256):
        assert gf256.add(a, a) == 0
 
 
def test_add_commutative():
    for a in range(0, 256, 13):
        for b in range(0, 256, 17):
            assert gf256.add(a, b) == gf256.add(b, a)
 
 
def test_sub_equals_add():
    for a in range(0, 256, 11):
        for b in range(0, 256, 13):
            assert gf256.sub(a, b) == gf256.add(a, b)
 
 
# ---------------------------------------------------------------------------
# Multiplication
# ---------------------------------------------------------------------------
 
def test_multiplicative_identity():
    for a in range(256):
        assert gf256.mul(a, 1) == a
        assert gf256.mul(1, a) == a
 
 
def test_zero_annihilates():
    for a in range(256):
        assert gf256.mul(a, 0) == 0
        assert gf256.mul(0, a) == 0
 
 
def test_mul_commutative():
    for a in range(0, 256, 7):
        for b in range(0, 256, 11):
            assert gf256.mul(a, b) == gf256.mul(b, a)
 
 
def test_mul_associative():
    samples = [0, 1, 2, 17, 42, 99, 128, 200, 255]
    for a in samples:
        for b in samples:
            for c in samples:
                lhs = gf256.mul(gf256.mul(a, b), c)
                rhs = gf256.mul(a, gf256.mul(b, c))
                assert lhs == rhs
 
 
def test_distributive():
    samples = [0, 1, 2, 17, 42, 99, 128, 200, 255]
    for a in samples:
        for b in samples:
            for c in samples:
                lhs = gf256.mul(a, gf256.add(b, c))
                rhs = gf256.add(gf256.mul(a, b), gf256.mul(a, c))
                assert lhs == rhs
 
 
def test_specific_product():
    # alpha * alpha^7 = alpha^8 = 0x1D
    assert gf256.mul(2, 128) == 0x1D
 
 
# ---------------------------------------------------------------------------
# Inverse and division
# ---------------------------------------------------------------------------
 
def test_inv_property():
    for a in range(1, 256):
        assert gf256.mul(a, gf256.inv(a)) == 1
 
 
def test_inv_involution():
    for a in range(1, 256):
        assert gf256.inv(gf256.inv(a)) == a
 
 
def test_inv_of_zero_raises():
    try:
        gf256.inv(0)
    except ZeroDivisionError:
        return
    raise AssertionError("inv(0) should have raised ZeroDivisionError")
 
 
def test_div_zero_numerator():
    for b in range(1, 256):
        assert gf256.div(0, b) == 0
 
 
def test_div_mul_roundtrip():
    for a in range(0, 256, 5):
        for b in range(1, 256, 7):
            assert gf256.mul(gf256.div(a, b), b) == a
 
 
def test_div_by_zero_raises():
    try:
        gf256.div(42, 0)
    except ZeroDivisionError:
        return
    raise AssertionError("div(_, 0) should have raised ZeroDivisionError")
 
 
# ---------------------------------------------------------------------------
# Power
# ---------------------------------------------------------------------------
 
def test_exp_zero():
    for a in range(1, 256):
        assert gf256.power(a, 0) == 1
 
 
def test_exp_one():
    for a in range(256):
        assert gf256.power(a, 1) == a
 
 
def test_fermat():
    # a^255 = 1 for all nonzero a
    for a in range(1, 256):
        assert gf256.power(a, 255) == 1
 
 
def test_squared_equals_mul():
    for a in range(256):
        assert gf256.power(a, 2) == gf256.mul(a, a)
 
 
def test_negative_is_inverse():
    for a in range(1, 256, 3):
        assert gf256.power(a, -1) == gf256.inv(a)
        assert gf256.power(a, -2) == gf256.inv(gf256.mul(a, a))
 
 
def test_zero_pow_zero():
    assert gf256.power(0, 0) == 1
 
 
def test_zero_pow_positive():
    assert gf256.power(0, 7) == 0
 
 
def test_zero_pow_negative_raises():
    try:
        gf256.power(0, -1)
    except ZeroDivisionError:
        return
    raise AssertionError("power(0, -1) should have raised ZeroDivisionError")
 
 
# ---------------------------------------------------------------------------
# Polynomial addition
# ---------------------------------------------------------------------------
 
def test_poly_add_basic():
    # (3 + 2x) + (1 + x^2)  =  2 + 2x + x^2  in GF(2^8)
    assert gf256.poly_add([3, 2], [1, 0, 1]) == [2, 2, 1]
 
 
def test_poly_add_cancellation_collapses_to_zero():
    assert gf256.poly_add([1, 2, 3], [1, 2, 3]) == [0]
 
 
def test_poly_add_identity():
    assert gf256.poly_add([1, 2, 3], [0]) == [1, 2, 3]
    assert gf256.poly_add([0], [1, 2, 3]) == [1, 2, 3]
 
 
def test_poly_add_strips_trailing_zeros():
    # (1 + x) + (2 + x) = 3 + 0  →  [3]
    assert gf256.poly_add([1, 1], [2, 1]) == [3]
 
 
# ---------------------------------------------------------------------------
# Polynomial scale
# ---------------------------------------------------------------------------
 
def test_poly_scale_by_zero():
    assert gf256.poly_scale([1, 2, 3], 0) == [0]
 
 
def test_poly_scale_by_one():
    assert gf256.poly_scale([1, 2, 3], 1) == [1, 2, 3]
 
 
def test_poly_scale_matches_poly_mul():
    for c in [1, 2, 42, 128, 255]:
        assert gf256.poly_scale([1, 2, 3, 4], c) == gf256.poly_mul([1, 2, 3, 4], [c])
 
 
# ---------------------------------------------------------------------------
# Polynomial multiplication
# ---------------------------------------------------------------------------
 
def test_poly_mul_by_zero_poly():
    assert gf256.poly_mul([1, 2, 3], [0]) == [0]
 
 
def test_poly_mul_by_one_poly():
    assert gf256.poly_mul([1, 2, 3], [1]) == [1, 2, 3]
 
 
def test_poly_mul_square_of_linear():
    # (1 + x)(1 + x) = 1 + 2x + x^2 = 1 + x^2 in char 2
    assert gf256.poly_mul([1, 1], [1, 1]) == [1, 0, 1]
 
 
def test_poly_mul_commutative():
    p = [1, 17, 42, 200]
    q = [3, 9, 128]
    assert gf256.poly_mul(p, q) == gf256.poly_mul(q, p)
 
 
def test_poly_mul_distributive():
    p = [1, 2, 3]
    q = [4, 5]
    r = [6, 7, 8]
    lhs = gf256.poly_mul(p, gf256.poly_add(q, r))
    rhs = gf256.poly_add(gf256.poly_mul(p, q), gf256.poly_mul(p, r))
    assert lhs == rhs
 
 
# ---------------------------------------------------------------------------
# Polynomial evaluation
# ---------------------------------------------------------------------------
 
def test_poly_eval_at_zero_returns_constant():
    assert gf256.poly_eval([42, 1, 2, 3], 0) == 42
 
 
def test_poly_eval_constant_poly():
    assert gf256.poly_eval([42], 17) == 42
 
 
def test_poly_eval_linear():
    # p(x) = 3 + x, evaluated at 5  →  3 XOR 5 = 6
    assert gf256.poly_eval([3, 1], 5) == 6
 
 
def test_poly_eval_root_of_linear_factor():
    # Polynomial (x + alpha) has root alpha (since -alpha = alpha in char 2).
    alpha = 2
    assert gf256.poly_eval([alpha, 1], alpha) == 0
 
 
def test_poly_eval_matches_manual():
    # p(x) = 1 + 2x + 3x^2 at x = alpha = 2.
    x = 2
    expected = 1 ^ gf256.mul(2, x) ^ gf256.mul(3, gf256.mul(x, x))
    assert gf256.poly_eval([1, 2, 3], x) == expected
 
 
# ---------------------------------------------------------------------------
# Polynomial divmod
# ---------------------------------------------------------------------------
 
def test_poly_divmod_exact_division_roundtrip():
    # Build a product, divide by one factor, expect zero remainder.
    a = [1, 2, 3]
    b = [4, 5]
    product = gf256.poly_mul(a, b)
    q, r = gf256.poly_divmod(product, b)
    assert q == a
    assert r == [0]
 
 
def test_poly_divmod_known_example():
    # (x^2 + 2x + 3) / (x + 1) in GF(2^8): quotient [3, 1], remainder [0].
    q, r = gf256.poly_divmod([3, 2, 1], [1, 1])
    assert q == [3, 1]
    assert r == [0]
 
 
def test_poly_divmod_general_roundtrip():
    # num == q*den + r must hold for arbitrary num, den.
    num = [7, 3, 1, 9, 42, 17]
    den = [2, 5, 1]
    q, r = gf256.poly_divmod(num, den)
    reconstructed = gf256.poly_add(gf256.poly_mul(q, den), r)
    assert reconstructed == gf256.poly_add(num, [0])
    assert len(r) < len(den)
 
 
def test_poly_divmod_div_by_zero_poly_raises():
    try:
        gf256.poly_divmod([1, 2, 3], [0])
    except ZeroDivisionError:
        return
    raise AssertionError("poly_divmod by [0] should have raised ZeroDivisionError")
 
 
def test_poly_divmod_numerator_smaller_than_denominator():
    # deg(num) < deg(den) → q = 0, r = num
    q, r = gf256.poly_divmod([1, 2], [1, 1, 1])
    assert q == [0]
    assert r == [1, 2]
 
 
# ---------------------------------------------------------------------------
# Polynomial derivative
# ---------------------------------------------------------------------------
 
def test_poly_deriv_constant_is_zero():
    assert gf256.poly_deriv([42]) == [0]
 
 
def test_poly_deriv_linear():
    # d/dx(3 + 7x) = 7
    assert gf256.poly_deriv([3, 7]) == [7]
 
 
def test_poly_deriv_x_squared_vanishes_in_char_2():
    # d/dx(x^2) = 0
    assert gf256.poly_deriv([0, 0, 1]) == [0]
 
 
def test_poly_deriv_x_cubed():
    # d/dx(x^3) = 3*x^2, and 3 mod 2 = 1, so result is x^2.
    assert gf256.poly_deriv([0, 0, 0, 1]) == [0, 0, 1]
 
 
def test_poly_deriv_mixed():
    # d/dx(5 + 3x + 7x^2 + 2x^3 + x^4):
    # odd-indexed coefficients (i=1, i=3) survive; others vanish.
    # → [3, 0, 2]
    assert gf256.poly_deriv([5, 3, 7, 2, 1]) == [3, 0, 2]
 
 
# ---------------------------------------------------------------------------
# Test harness
# ---------------------------------------------------------------------------
 
def _run_all_tests():
    tests = sorted(
        (name, fn)
        for name, fn in globals().items()
        if name.startswith("test_") and callable(fn)
    )
    passed = 0
    failed = []
    for name, fn in tests:
        try:
            fn()
        except Exception:
            failed.append((name, traceback.format_exc()))
            print(f"  FAIL  {name}")
        else:
            passed += 1
            print(f"  ok    {name}")
    print()
    print(f"{passed}/{len(tests)} tests passed")
    if failed:
        print()
        print("FAILURES:")
        for name, tb in failed:
            print(f"\n--- {name} ---")
            print(tb)
        return 1
    return 0
 
 
if __name__ == "__main__":
    sys.exit(_run_all_tests())