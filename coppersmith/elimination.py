from __future__ import annotations

from fractions import Fraction

# 消元与结果式工具（不依赖外部库）
# - 针对二元多项式的“按 y 视作一元”结果式 R(x)
# - 通过多点取值与插值恢复 R(x)（避免在 Z[x] 上直接行列式）

BivarFrac = dict[tuple[int, int], Fraction]  # (ix,iy) -> Fraction


# ----------------- 基础工具 -----------------


def deg_x_bivar_frac(F: BivarFrac) -> int:
    return max((ix for (ix, _iy) in F.keys()), default=-1)


def deg_y_bivar_frac(F: BivarFrac) -> int:
    return max((iy for (_ix, iy) in F.keys()), default=-1)


def bivar_frac_eval_x_get_univar_y(F: BivarFrac, x0: int) -> list[Fraction]:
    """将 F(x,y) 在 x=x0 处专化，得到关于 y 的一元多项式（升幂系数，Fraction）。"""
    dy = deg_y_bivar_frac(F)
    coeffs = [Fraction(0) for _ in range(max(dy, -1) + 1)]
    for (ix, iy), c in F.items():
        coeffs[iy] += c * (x0**ix)
    # 去掉尾部 0
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def clear_denominators(fr_coeffs: list[Fraction]) -> tuple[list[int], int]:
    """将分数系数清分母，返回整数系数与公共分母 D（使得 int_coeffs = fr_coeffs * D）。"""
    from math import gcd

    if not fr_coeffs:
        return [], 1

    # 计算最小公倍数作为公共分母
    def lcm(a: int, b: int) -> int:
        return (a * b) // gcd(a, b)

    # 依次累乘 LCM，避免 reduce 的类型歧义
    D = 1
    for frac in fr_coeffs:
        D = lcm(D, frac.denominator)
    int_coeffs = [int(a * D) for a in fr_coeffs]

    # 规范化去公因子（迭代 gcd）
    g = 0
    for v in int_coeffs:
        g = gcd(g, abs(v))
    if g > 1:
        int_coeffs = [v // g for v in int_coeffs]
        D //= g
    return int_coeffs, D


# ----------------- 整数行列式（Bareiss）与结果式 -----------------


def det_bareiss_int(A: list[list[int]]) -> int:
    """Bareiss 无分式消元法，整数矩阵行列式。A 不会被修改。"""
    n = len(A)
    if n == 0:
        return 1
    M = [row[:] for row in A]
    sign = 1
    prev = 1
    for k in range(n - 1):
        # 选主元（若为 0，尝试行交换）
        if M[k][k] == 0:
            pivot_row = None
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    pivot_row = r
                    break
            if pivot_row is None:
                return 0
            M[k], M[pivot_row] = M[pivot_row], M[k]
            sign = -sign
        pivot = M[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                # M[i][j] = (M[i][j]*pivot - M[i][k]*M[k][j]) / prev
                M[i][j] = (M[i][j] * pivot - M[i][k] * M[k][j]) // prev
            M[i][k] = 0
        prev = pivot
    return sign * M[n - 1][n - 1]


def sylvester_matrix_int(p: list[int], q: list[int]) -> list[list[int]]:
    """关于 y 的 Sylvester 矩阵（整数），p,q 为升幂系数。"""
    dp = len(p) - 1
    dq = len(q) - 1
    n = dp + dq
    # 行 n，列 n
    S = [[0 for _ in range(n)] for _ in range(n)]
    # 前 dq 行放 p 的移位
    for row in range(dq):
        for i, ci in enumerate(p):
            col = row + i
            if col < n:
                S[row][col] = ci
    # 后 dp 行放 q 的移位
    for row in range(dp):
        r = dq + row
        for i, ci in enumerate(q):
            col = row + i
            if col < n:
                S[r][col] = ci
    return S


def resultant_int_y(p: list[int], q: list[int]) -> int:
    if not p or not q:
        return 0
    if all(v == 0 for v in p) or all(v == 0 for v in q):
        return 0
    S = sylvester_matrix_int(p, q)
    return det_bareiss_int(S)


# ----------------- 插值重建 R(x) -----------------


def lagrange_interpolate(points: list[tuple[int, int]]) -> list[int]:
    """拉格朗日插值（整数输出）。输入点 x 均不相同。返回升幂整数系数。"""
    # 使用 Fraction 做精确计算
    xs = [Fraction(x) for x, _ in points]
    ys = [Fraction(y) for _, y in points]
    # 多项式累加法
    coeffs = [Fraction(0)]
    for i in range(len(points)):
        # 构造基函数 Li(x)
        num = [Fraction(1)]  # 多项式 1
        den = Fraction(1)
        xi = xs[i]
        for j in range(len(points)):
            if i == j:
                continue
            xj = xs[j]
            # num *= (x - xj)
            num = poly_mul_frac(num, [-xj, Fraction(1)])
            den *= xi - xj
        Li = [c / den for c in num]
        term = [c * ys[i] for c in Li]
        coeffs = poly_add_frac(coeffs, term)
    # 规约为整数系数
    return [int(c) for c in normalize_fraction_coeffs(coeffs)]


def poly_add_frac(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    n = max(len(a), len(b))
    out = [Fraction(0) for _ in range(n)]
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += b[i]
    return out


def poly_mul_frac(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0) for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def normalize_fraction_coeffs(a: list[Fraction]) -> list[Fraction]:
    # 清尾零
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    # 转整数
    from math import gcd

    D = 1
    for c in a:
        D = (D * c.denominator) // gcd(D, c.denominator)
    ints = [int(c * D) for c in a]
    g = 0
    for v in ints:
        g = gcd(g, abs(v))
    if g > 1:
        ints = [v // g for v in ints]
        D //= g
    return [Fraction(v, 1) for v in ints]


def eval_int_poly(coeffs: list[int], x: int) -> int:
    total = 0
    for i, ci in enumerate(coeffs):
        total += ci * (x**i)
    return total


def resultant_in_x_by_interpolation(G1: BivarFrac, G2: BivarFrac, X: int, Y: int) -> list[int]:
    """计算关于 y 的结果式 R(x)，返回整数多项式（升幂系数）。
    方法：选取一批 x0 点，专化得到一元多项式 P1(y),P2(y)（分数系数），清分母为整数后求整数结果式；
    收集 (x0, R'(x0)) 点，用拉格朗日插值重建 R'(x)；该 R'(x) 与真 R(x) 仅差一个非零常数因子。
    """
    # 度数上界：deg_x(R) ≤ deg_y(G1)*deg_x(G2) + deg_y(G2)*deg_x(G1)
    dxy1 = deg_x_bivar_frac(G1)
    dyy1 = deg_y_bivar_frac(G1)
    dxy2 = deg_x_bivar_frac(G2)
    dyy2 = deg_y_bivar_frac(G2)
    deg_bound = max(0, dyy1 * dxy2 + dyy2 * dxy1)

    # 选择采样点（尽量小整数，覆盖 [-X-2, X+2] 范围）
    samples: list[tuple[int, int]] = []
    needed = deg_bound + 1
    # 构造候选点列表
    cand = list(range(-(X + 3), X + 4))
    # 去重并过滤重复
    used = set()
    for x0 in cand:
        if x0 in used:
            continue
        used.add(x0)
        # 专化并清分母
        p1_y = bivar_frac_eval_x_get_univar_y(G1, x0)
        p2_y = bivar_frac_eval_x_get_univar_y(G2, x0)
        int_p1, D1 = clear_denominators(p1_y)
        int_p2, D2 = clear_denominators(p2_y)
        if len(int_p1) == 0 or len(int_p2) == 0:
            continue
        # 求整数结果式
        r_val = resultant_int_y(int_p1, int_p2)
        samples.append((x0, r_val))
        if len(samples) >= needed:
            break
    if len(samples) < needed:
        # 退化：直接返回常数多项式
        return [0, 0, 1][:1]  # [0]

    coeffs = lagrange_interpolate(samples)
    # 归一化：去除系数的最大公因子
    from math import gcd

    g = 0
    for v in coeffs:
        g = gcd(g, abs(v))
    if g > 1:
        coeffs = [v // g for v in coeffs]
    return coeffs
