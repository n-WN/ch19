from __future__ import annotations

from fractions import Fraction

# 简单整数 LLL 实现（列向量基或行向量基的一致性）
# 这里使用“行向量”为基，输入为矩阵 rows: List[List[int]]
# 输出同维度的约化基（行向量）

Vector = list[Fraction]
Matrix = list[Vector]


def dot(a: Vector, b: Vector) -> Fraction:
    # 指定 Fraction(0) 作为起始值，确保返回类型为 Fraction
    # 注意：ty 对 zip(strict=...) 的类型支持较保守，这里避免使用 strict 以通过类型检查
    total = Fraction(0)
    for i in range(len(a)):
        total += a[i] * b[i]
    return total


def gram_schmidt(B: Matrix) -> tuple[Matrix, list[Fraction]]:
    """Compute Gram–Schmidt orthogonalization of rows in B.

    Returns:
      (U, norms2): U is the orthogonalized basis; norms2 are squared norms of U[i].
    """
    n = len(B)
    m = len(B[0]) if n else 0
    U: Matrix = [[Fraction(0) for _ in range(m)] for _ in range(n)]
    mu: list[list[Fraction]] = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    Bstar_norm2: list[Fraction] = [Fraction(0) for _ in range(n)]

    for i in range(n):
        U[i] = list(B[i])
        for j in range(i):
            mu[i][j] = dot(B[i], U[j]) / Bstar_norm2[j] if Bstar_norm2[j] != 0 else Fraction(0)
            # U[i] -= mu[i][j] * U[j]
            for k in range(m):
                U[i][k] -= mu[i][j] * U[j][k]
        Bstar_norm2[i] = dot(U[i], U[i])
    return U, [Bstar_norm2[i] for i in range(n)]


def size_reduce(B: Matrix, mu: list[list[Fraction]], k: int, ell: int) -> None:
    r = int(round(mu[k][ell]))
    if r != 0:
        for j in range(len(B[k])):
            B[k][j] -= r * B[ell][j]
        mu[k][ell] -= r
        for j in range(ell):
            mu[k][j] -= r * mu[ell][j]


def lll_reduction(B_int: list[list[int]], delta: Fraction = Fraction(3, 4)) -> list[list[int]]:
    # 转为 Fraction
    B: Matrix = [[Fraction(x) for x in row] for row in B_int]
    n = len(B)
    if n == 0:
        return []
    _m = len(B[0])
    k = 1
    U, Bstar_norm2 = gram_schmidt(B)
    mu: list[list[Fraction]] = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i):
            mu[i][j] = dot(B[i], U[j]) / Bstar_norm2[j] if Bstar_norm2[j] != 0 else Fraction(0)

    while k < n:
        # size reduction 对 j = k-1
        size_reduce(B, mu, k, k - 1)
        U, Bstar_norm2 = gram_schmidt(B)
        for i in range(n):
            for j in range(i):
                mu[i][j] = dot(B[i], U[j]) / Bstar_norm2[j] if Bstar_norm2[j] != 0 else Fraction(0)

        # Lovasz 条件
        left = Bstar_norm2[k]
        right = (delta - mu[k][k - 1] ** 2) * Bstar_norm2[k - 1]
        if left >= right:
            # 对 j = k-2, ..., 0 做 size reduction
            for j in reversed(range(k - 1)):
                size_reduce(B, mu, k, j)
            k += 1
        else:
            # 交换 b_k 和 b_{k-1}
            B[k], B[k - 1] = B[k - 1], B[k]
            U, Bstar_norm2 = gram_schmidt(B)
            for i in range(n):
                for j in range(i):
                    mu[i][j] = (
                        dot(B[i], U[j]) / Bstar_norm2[j] if Bstar_norm2[j] != 0 else Fraction(0)
                    )
            k = max(k - 1, 1)

    # 转回 int（四舍五入）
    return [[int(x) for x in row] for row in B]
