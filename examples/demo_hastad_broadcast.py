#!/usr/bin/env python3
from __future__ import annotations

import random

# Hastad 广播攻击（e=3）：同一消息 m 在不同互素模数 Ni 上加密 ci = m^3 mod Ni
# 在无填充情况下，若 m^3 < N1*N2*N3，则 CRT 合并得到 C ≡ m^3 (mod N123)
# 且 0 ≤ C < N123，于是 m = ⌊C^{1/3}⌉。这里完整实现：手写 CRT 与整数立方根。


def egcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = egcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)


def inv_mod(a: int, m: int) -> int:
    g, x, y = egcd(a, m)
    if g != 1:
        raise ValueError("not invertible")
    return x % m


def crt(residues: list[int], moduli: list[int]) -> tuple[int, int]:
    # 返回 (r, M)，满足 r ≡ residues[i] (mod moduli[i])，M = Π moduli
    M = 1
    for n in moduli:
        M *= n
    r = 0
    # 避免 strict 以兼容 ty 的类型检查
    for idx, ai in enumerate(residues):
        ni = moduli[idx]
        Mi = M // ni
        ti = inv_mod(Mi % ni, ni)
        r = (r + ai * Mi * ti) % M
    return r, M


def integer_cuberoot_floor(n: int) -> int:
    # 整数立方根下取整
    if n < 0:
        raise ValueError("n must be >= 0")
    x = int(round(n ** (1 / 3)))
    # 校正
    while (x + 1) ** 3 <= n:
        x += 1
    while x**3 > n:
        x -= 1
    return x


def main() -> None:
    random.seed(999)
    e = 3

    # 生成三个互素模数（小规模演示）
    def gen_prime(bits: int) -> int:
        from math import isqrt

        def is_prime(n: int) -> bool:
            if n < 2:
                return False
            if n % 2 == 0:
                return n == 2
            r = isqrt(n)
            f = 3
            while f <= r:
                if n % f == 0:
                    return False
                f += 2
            return True

        while True:
            x = random.getrandbits(bits) | 1
            if is_prime(x):
                return x

    N1 = gen_prime(20) * gen_prime(20)
    N2 = gen_prime(20) * gen_prime(20)
    N3 = gen_prime(20) * gen_prime(20)

    # 选择小消息 m，使 m^3 < N1*N2*N3
    Mprod = N1 * N2 * N3
    X = int(round(Mprod ** (1 / 3)))
    m_true = max(2, X // 4)

    c1 = pow(m_true, e, N1)
    c2 = pow(m_true, e, N2)
    c3 = pow(m_true, e, N3)

    C, M = crt([c1, c2, c3], [N1, N2, N3])
    # 此时 C ≡ m^3 (mod M) 且 m^3 < M，故 C == m^3
    m_rec = integer_cuberoot_floor(C)

    print(
        {
            "N1": N1,
            "N2": N2,
            "N3": N3,
            "m_true": m_true,
            "C": C,
            "m_rec": m_rec,
        }
    )

    assert m_rec == m_true, "Hastad 广播攻击失败；可增大位数重新生成"
    print("OK: Hastad 广播攻击复现成功")


if __name__ == "__main__":
    main()
