import math

"""
# Gauss-Legendre nodes and weights for n = 32 on [-1, 1].
# Source: https://pomax.github.io/bezierinfo/legendre-gauss.html
# Listed as (node, weight) pairs covering only the positive half of the symmetric distribution.
# Organized by generative AI
"""
_GL32 = [
    (0.0483076656877383, 0.0965400885147278),
    (0.1444719615827965, 0.0956387200792749),
    (0.2392873622521371, 0.0938443990808046),
    (0.3318686022821277, 0.0911738786957639),
    (0.4213512761306353, 0.0876520930044038),
    (0.5068999089322294, 0.0833119242269467),
    (0.5877157572407623, 0.0781938957870703),
    (0.6630442669302152, 0.0723457941088485),
    (0.7321821187402897, 0.0658222227763618),
    (0.7944837959679424, 0.0586840934785355),
    (0.8493676137325700, 0.0509980592623762),
    (0.8963211557660521, 0.0428358980222267),
    (0.9349060759377397, 0.0342738629130214),
    (0.9647622555875064, 0.0253920653092621),
    (0.9856115115452684, 0.0162743947309057),
    (0.9972638618494816, 0.0070186100094701),
]

def gauss_legendre(f, a, b, n=32):
    # Numerically integrate f over [a, b] using n-point Gauss-Legendre
    half_width = (b - a) / 2.0
    midpoint = (b + a) / 2.0

    total = 0.0
    for node, weight in _GL32:
        # Each entry in _GL32 represents two symmetric nodes (+/- node). Evaluate f at both and accumulate
        t_pos = half_width * node + midpoint
        t_neg = half_width * (-node) + midpoint
        total += weight * (f(t_pos) + f(t_neg))

    return half_width * total

_INV_SQRT_2PI = 1.0 / math.sqrt(2.0 * math.pi)

def q_function(x):
    """
    Compute Q(x) = P(Z > x) for standard normal Z, via numerical integration

    For x >= 0, integrates exp(-u^2/2) from x to x + 20 and scales by 1/sqrt(2*pi). For x < 0, uses the symmetry Q(-x) = 1 - Q(x).

    Used to verify the closed-form uncoded BPSK BER expression BER = Q(sqrt(2 * Eb/N0)) against numerical integration, satisfying proposal R4
    """
    if x < 0:
        return 1.0 - q_function(-x)

    integrand = lambda u: math.exp(-0.5 * u * u)
    integral = gauss_legendre(integrand, x, x + 20.0, n=32)
    return _INV_SQRT_2PI * integral

def central_diff(values, h):
    """
    Approximate the derivative of a sampled function

    For a function with uniformly spaced  points at spacing h, the central difference is 
        f'(x_i) ~= (y_{i+1} - y_{i-1}) / (2h)
    At the endpoints i = 0 and i = n-1 the central formula is undefined (it would index outside the array), so we fall back to one-sided forward/backward differences
        f'(x_0) ~= (y_1 - y_0) / h (forward)
        f'(x_{n-1}) ~= (y_{n-1} - y_{n-2}) / h (backward)
    """
    n = len(values)
    derivatives = [0.0] * n

    # Forward difference at the left boundary.
    derivatives[0] = (values[1] - values[0]) / h

    # Central differences at interior points.
    for i in range(1, n - 1):
        derivatives[i] = (values[i + 1] - values[i - 1]) / (2.0 * h)

    # Backward difference at the right boundary.
    derivatives[-1] = (values[-1] - values[-2]) / h

    return derivatives