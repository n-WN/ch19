from __future__ import annotations

# 多项式用 dict[int, int] 存储：{幂次: 系数}，系数为 int，自动规范化（去零）

Poly = dict[int, int]


def normalize(p: Poly) -> Poly:
    """Return a new polynomial with zero coefficients removed."""
    return {k: v for k, v in p.items() if v != 0}


def from_coeffs(coeffs: list[int]) -> Poly:
    """Create a polynomial from ascending coefficients.

    Args:
      coeffs: [c0, c1, c2, ...] meaning c0 + c1*x + c2*x^2 + ...
    """
    p: Poly = {}
    for i, c in enumerate(coeffs):
        if c:
            p[i] = int(c)
    return p


def to_coeffs(p: Poly) -> list[int]:
    """Convert dict-based poly to dense ascending coefficients list."""
    if not p:
        return [0]
    deg = max(p)
    out = [0] * (deg + 1)
    for k, v in p.items():
        out[k] = v
    return out


def degree(p: Poly) -> int:
    """Return the degree of the polynomial, or -1 for the zero polynomial."""
    return max(p.keys(), default=-1)


def add(a: Poly, b: Poly) -> Poly:
    """Return a + b."""
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
        if out[k] == 0:
            del out[k]
    return out


def sub(a: Poly, b: Poly) -> Poly:
    """Return a - b."""
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) - v
        if out[k] == 0:
            del out[k]
    return out


def mul(a: Poly, b: Poly) -> Poly:
    """Return a * b (convolution)."""
    out: Poly = {}
    for i, ai in a.items():
        for j, bj in b.items():
            out[i + j] = out.get(i + j, 0) + ai * bj
            if out[i + j] == 0:
                del out[i + j]
    return out


def mul_xk(a: Poly, k: int) -> Poly:
    """Return a * x^k.

    Raises:
      ValueError: if k < 0
    """
    if k < 0:
        raise ValueError("k must be >= 0")
    out: Poly = {}
    for i, ai in a.items():
        out[i + k] = ai
    return out


def scale(a: Poly, c: int) -> Poly:
    """Return c * a."""
    if c == 0:
        return {}
    return {i: ai * c for i, ai in a.items()}


def pow_poly(a: Poly, e: int) -> Poly:
    """Return a**e using binary exponentiation.

    Raises:
      ValueError: if e < 0
    """
    if e < 0:
        raise ValueError("exponent must be >= 0")
    out: Poly = {0: 1}
    base = dict(a)
    ee = e
    while ee > 0:
        if ee & 1:
            out = mul(out, base)
        base = mul(base, base)
        ee >>= 1
    return out


def eval_at(a: Poly, x: int) -> int:
    """Evaluate a(x).

    Note: Horner's scheme may be faster; here we favor clarity.
    """
    # Horner 法也可，这里直接用幂展开
    total = 0
    for i, ai in a.items():
        total += ai * (x**i)
    return total


def mod_poly(a: Poly, m: int) -> Poly:
    """Reduce coefficients modulo m (m>0).

    Raises:
      ValueError: if m <= 0
    """
    if m <= 0:
        raise ValueError("modulus must be positive")
    out: Poly = {}
    for i, ai in a.items():
        r = ai % m
        if r:
            out[i] = r
    return out
