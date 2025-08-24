#!/usr/bin/env python3
from __future__ import annotations

import random
from math import isqrt

import pytest

from coppersmith.bivar import Bivar
from coppersmith.bivariate import try_find_small_roots_bivar

# 基础正确性与性能回归测试
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


@pytest.mark.parametrize("seed,bits,X", [(1, 14, 256), (2, 14, 256)])
def test_univariate_small_root(seed: int, bits: int, X: int) -> None:
    random.seed(seed)
    p = gen_prime(bits)
    q = gen_prime(bits)
    N = p * q
    r = random.randrange(-X // 2, X // 2)
    c = (-(r * r)) % N
    f_coeffs = [c, 0, 1]
    roots = find_small_roots_univariate(f_coeffs, N=N, X=X, m=3, t=3)
    print({"case": "univar", "N_bits": bits * 2, "X": X, "r": r, "roots": roots})
    assert r in roots


def test_bivariate_constructed_root() -> None:
    random.seed(7)
    p = 499
    q = 547
    N = p * q
    X = 24
    Y = 24
    r = -2
    s = -8
    c = (-(r * r + s)) % N
    F: Bivar = {(2, 0): 1, (0, 1): 1, (0, 0): c}
    roots = try_find_small_roots_bivar(F, N=N, X=X, Y=Y, m=2, tx=2, ty=2)
    print({"case": "bivar", "N": N, "X": X, "Y": Y, "true": (r, s), "roots": roots[:10]})
    assert (r, s) in roots
