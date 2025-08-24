from __future__ import annotations

from fractions import Fraction

from .lll import lll_reduction
from .poly import Poly, degree, from_coeffs, mul_xk, pow_poly, scale

# 教学版：单变量 Coppersmith 小根方法（基础版，Howgrave-Graham 变体）
# 输入：
# - f_coeffs: 升幂系数列表（整数），表示 f(x) ∈ Z[x]
# - N: 模数（正整数）
# - X: 对小根的界，期望 |x0| < X
# - m, t: 构造格的参数（常见取 m≈deg(f)，t≈deg(f) 或更小）
# 输出：
# - 候选整数根列表（去重），每个 r 满足 |r|<X 且 f(r)≡0 (mod N)


def construct_lattice(
    f_coeffs: list[int], N: int, X: int, m: int, t: int
) -> tuple[list[list[int]], int]:
    """构造格基：
    - 对 i = 0..m-1, j = 0..d-1:  N^{m-i} * x^j * f(x)^i
    - 对 i = 0..t-1:              x^i * f(x)^m
    然后进行列缩放：列 k 乘以 X^k，相当于对变量替换 x -> X·x
    返回：整数矩阵 B 以及列数 ncols
    """
    f = from_coeffs(f_coeffs)
    d = degree(f)

    polys: list[Poly] = []
    # i = 0..m-1 层
    for i in range(m):
        fi = pow_poly(f, i)
        scale_N = pow(N, m - i)
        fi_scaled = scale(fi, scale_N)
        polys.extend(mul_xk(fi_scaled, j) for j in range(d))
    # f(x)^m 的 t 个移位
    f_m = pow_poly(f, m)
    polys.extend(mul_xk(f_m, i) for i in range(t))

    # 列缩放并转为整数矩阵
    max_deg = max((degree(p) for p in polys), default=0)
    ncols = max_deg + 1
    B: list[list[int]] = []
    for p in polys:
        row = [0] * ncols
        for k, ak in p.items():
            row[k] = ak * pow(X, k)
        B.append(row)
    return B, ncols


def eval_unscaled_row_at(row: list[int], r: int, X: int) -> Fraction:
    """Evaluate scaled polynomial row at integer r after unscaling by powers of X."""
    """行向量 row 是缩放后多项式（替换 x->X·x）的系数。
    反缩放在点 r 处的值：sum (row[k]/X^k) * r^k = 以原变量评估。
    返回 Fraction 以避免精度问题。
    """
    total = Fraction(0)
    # 缓存 X 和 r 的幂
    X_pow = 1
    r_pow = 1
    for _k, ck in enumerate(row):
        if ck:
            total += Fraction(ck, X_pow) * r_pow
        X_pow *= X
        r_pow *= r
    return total


def find_small_roots_univariate(
    f_coeffs: list[int], N: int, X: int, m: int = 3, t: int = 3
) -> list[int]:
    """Find small roots |r|<X of f(x) ≡ 0 (mod N) using lattice/LLL.

    Args:
      f_coeffs: ascending integer coefficients of f
      N: modulus (>0)
      X: search bound (>0)
      m,t: lattice parameters
    Returns:
      Sorted list of integer roots r with |r|<X and f(r)≡0 (mod N)
    """
    if N <= 0:
        raise ValueError("N must be positive")
    if X <= 0:
        return []
    if not f_coeffs:
        return []

    B, _ = construct_lattice(f_coeffs, N, X, m, t)
    Bref = lll_reduction(B)

    candidates = set()
    # 取前若干短向量尝试
    for row in Bref[: min(len(Bref), 12)]:
        for r in range(-X + 1, X):
            val = eval_unscaled_row_at(row, r, X)
            if val == 0:
                # 验证 f(r) ≡ 0 (mod N)
                fr = 0
                for i, ai in enumerate(f_coeffs):
                    fr += ai * pow(r, i)
                if fr % N == 0:
                    candidates.add(r)
    return sorted(candidates)
