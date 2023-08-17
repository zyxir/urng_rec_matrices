"""
Search and store matrices for uniform random number generation.

Based on \"High Quality Uniform Random Number Generation Through LUT Optimized
Linear Recurrences\" by David B.  Thomas and Wayne Luk in 2005.
"""

import json
import time
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import sympy
from tqdm import tqdm


def gen_matrix(k: int, t: int, e: int) -> np.ndarray:
    """Generate an special `k`x`k` matrix.

    The matrix have exactly `t` ones in each row and column, with `e`
    exceptions: in the `e` rows there are only `t-1` ones.
    """
    result = None
    while result is None:
        result = _gen_matrix(k, t, e)
    return result


def _gen_matrix(k: int, t: int, e: int) -> Optional[np.ndarray]:
    """Try to generate one special matrix."""
    col_sums = np.zeros(k)
    avail_cols = list(range(k))
    exceptions = np.random.choice(k, e)
    mat = np.zeros((k, k), dtype=int)
    try:
        for i in range(k):
            t_actual = t - 1 if i in exceptions else t
            indices = np.random.choice(avail_cols, t_actual, replace=False)
            mat[i, indices] = 1
            for j in indices:
                col_sums[j] += 1
                if col_sums[j] == t:
                    avail_cols.remove(j)
        return mat
    except ValueError:
        return None


def search_matrix(k: int, t: int) -> Tuple[np.ndarray, Dict[str, Union[int, float]]]:
    """Search for a matrix that satisfies these specific criteria:

        1. Has a shape of `k`x`k`.

        2. Each row and column have `t` ones, except 1 of each can have `t-1`
           ones.  Other elements are zeros.

        3. Have a primitive characteristic polynomial.

    Return the found matrix with additonal information.
    """
    if t >= k:
        raise ValueError("t must be smaller than k!")

    attempts = 0
    det_rejects = 0
    irre_rejects = 0
    prim_rejects = 0
    max_attempts = 10000

    start_time = time.time()
    while attempts < max_attempts:
        attempts = attempts + 1
        # Generate new random matrix.
        mat = gen_matrix(k, t, 1)
        # Check determinant.
        det = np.linalg.det(mat)
        if det == 0:
            det_rejects += 1
            continue
        # Check CP reducibility.
        x = sympy.Symbol("x")
        cp = sympy.Poly(sympy.Matrix(mat).charpoly(x).as_expr(), domain=sympy.QQ)
        if not cp.is_irreducible:
            irre_rejects += 1
            continue
        # Check CP primitivity.
        if cp.is_primitive:
            break
        else:
            prim_rejects += 1
    end_time = time.time()
    elapsed_time = end_time - start_time

    info = {
        "attempts": attempts,
        "det_rejects": det_rejects,
        "irre_rejects": irre_rejects,
        "prim_rejects": prim_rejects,
        "elapsed_time": elapsed_time,
    }
    return (mat, info)


def store_matrix(mat: np.ndarray, k: int) -> List[List[int]]:
    """Store matrix in a sparse format.

    The `k`x`k` matrix will be stored a list of size `k`.  The `i`th element is
    a list of the indices of 1s in the `i`th row.
    """
    result = []
    for i in range(k):
        elem = []
        for j in range(k):
            if mat[i, j] == 1:
                elem.append(j)
        result.append(elem)
    return result


if __name__ == "__main__":
    # The list of `k`s and `t`s (numbers of inputs of LUT) to calculate.
    ks = range(8, 256 + 8, 8)
    ts = [4, 6, 8]

    # Assume a `k**2` time complexity, estimate the total progress.
    total = (np.asarray(ks)**2 * len(ts)).sum()
    progress_bar = tqdm(total=total, desc="Processing")

    # Generate all matrices.
    data = dict()
    for k in ks:
        elem = dict()
        for t in ts:
            if t < k:
                mat, _ = search_matrix(k, t)
                elem[t] = store_matrix(mat, k)
            progress_bar.update(k**2)
        data[k] = elem

    # Save the matrices as a JSON file.
    fname = "urng_rec_matrices.json"
    with open(fname, "w") as fp:
        content = json.dumps(data, indent=4)
        fp.write(content)
        print(f"{fname} successfully written.")
