#!/usr/bin/env python3
from __future__ import annotations

import random
import time
from math import isqrt

import pytest

from coppersmith.univariate import find_small_roots_univariate


def gen_prime(bits: int) -> int:
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


@pytest.mark.parametrize("bits,X", [(16, 256), (16, 384)])
def test_univariate_performance(bits: int, X: int) -> None:
    random.seed(42)
    p = gen_prime(bits)
    q = gen_prime(bits)
    N = p * q
    r = random.randrange(-X // 2, X // 2)
    c = (-(r * r)) % N
    f_coeffs = [c, 0, 1]

    t0 = time.perf_counter()
    roots = find_small_roots_univariate(f_coeffs, N=N, X=X, m=3, t=3)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    print(
        {"case": "perf_univar", "X": X, "elapsed_sec": round(elapsed, 4), "roots_len": len(roots)}
    )
    # 正确性
    assert r in roots
    # 性能阈值（教学参数下 0.4s 内）
    assert elapsed < 0.4
