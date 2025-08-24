from __future__ import annotations

# 二元多项式：map[(ix, iy)] -> int
Bivar = dict[tuple[int, int], int]


def normalize(p: Bivar) -> Bivar:
    return {k: v for k, v in p.items() if v != 0}


def add(a: Bivar, b: Bivar) -> Bivar:
    out = dict(a)
    for (i, j), v in b.items():
        out[(i, j)] = out.get((i, j), 0) + v
        if out[(i, j)] == 0:
            del out[(i, j)]
    return out


def scale(a: Bivar, c: int) -> Bivar:
    if c == 0:
        return {}
    return {(i, j): v * c for (i, j), v in a.items()}


def mul(a: Bivar, b: Bivar) -> Bivar:
    out: Bivar = {}
    for (i1, j1), v1 in a.items():
        for (i2, j2), v2 in b.items():
            key = (i1 + i2, j1 + j2)
            out[key] = out.get(key, 0) + v1 * v2
            if out[key] == 0:
                del out[key]
    return out


def pow_bivar(a: Bivar, e: int) -> Bivar:
    if e < 0:
        raise ValueError("exponent must be >= 0")
    out: Bivar = {(0, 0): 1}
    base = dict(a)
    ee = e
    while ee > 0:
        if ee & 1:
            out = mul(out, base)
        base = mul(base, base)
        ee >>= 1
    return out


def shift_x(a: Bivar, k: int) -> Bivar:
    if k < 0:
        raise ValueError("k must be >= 0")
    return {(i + k, j): v for (i, j), v in a.items()}


def shift_y(a: Bivar, k: int) -> Bivar:
    if k < 0:
        raise ValueError("k must be >= 0")
    return {(i, j + k): v for (i, j), v in a.items()}


def degree_x(a: Bivar) -> int:
    return max((i for (i, _j) in a.keys()), default=-1)


def degree_y(a: Bivar) -> int:
    return max((j for (_i, j) in a.keys()), default=-1)


def eval_at(a: Bivar, x: int, y: int) -> int:
    total = 0
    for (i, j), v in a.items():
        total += v * (x**i) * (y**j)
    return total
