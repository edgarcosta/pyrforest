from .rforest import remainder_forest

from sage.matrix.constructor import Matrix
from sage.rings.fast_arith import prime_range
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

def batch_factorial(n, e, gamma):
    r"""
    Amortized computation of factorials.

    INPUT:

     - ``n``: bound on primes
     - ``e``: order of the series expansion
     - ``gamma``: a rational number

    OUTPUT:

    A dict whose value at a prime `p` equals `(\lceil \gamma*p\rceil - 1)! \pmod{p^e}`.

    EXAMPLES::

        sage: from pyrforest.batch import *
        sage: batch_factorial(10, 2, 1/2)
        {2: 1, 3: 1, 5: 2, 7: 6}
    """
    R = ZZ["k"]
    M = Matrix(R, 1, 1, [R.gen()])
    a, b = QQ(gamma).as_integer_ratio()
    k = lambda p, a=a, b=b: -(-a * p // b)
    m = lambda p, e=e: p**e

    ans = remainder_forest(M, m, k, kbase=1, indices=prime_range(n))
    return {p: ans[p][0, 0] for p in prime_range(n)}

def batch_harmonic(n, e, gamma, j, proj=False):
    r"""
    Amortized computation of harmonic sums.

    INPUT:

     - ``n``: bound on primes
     - ``e``: order of the series expansion
     - ``gamma``: a rational number
     - ``j``: a positive integer

    OUTPUT:

    A dict whose value at a prime `p` is the truncated harmonic sum
    `\sum_{k=1}^{\lceil \gamma*p \rceil - 1} k^{-j} \pmod{p^e}`.

    If ``proj`` is True, instead return pairs `(x,y)` representing `x/y \pmod{p^e}`.
    """
    # Represent summation as a 2x2 matrix.
    R = ZZ["k"]
    y = R.gen()
    M = Matrix(R, [[y**j, 0], [1, y**j]])

    a, b = QQ(gamma).as_integer_ratio()
    k = lambda p, a=a, b=b: -(-a * p // b)
    m = lambda p, e=e: p**e
    V = Matrix(ZZ, [[0, 1]])

    ans = remainder_forest(M, m, k, kbase=1, indices=prime_range(n), V=V)
    if proj:
        return ans
    return {p: mat[0, 0] * mat[0, 1].inverse_mod(p**e) for p, mat in ans.items()}
