from .rforest import remainder_forest

from sage.matrix.constructor import Matrix
from sage.rings.fast_arith import prime_range
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

def batch_factorial(n, e, gamma):
    r"""
    Amortized computation of factorials.

    INPUT::

     - `n`: bound on primes
     - `e`: order of the series expansion
     - `gamma`: a rational number

    OUTPUT::
       A dict whose value at a prime `p` equals (ceil(gamma*p)-1)! mod p^e.
    """
    R = ZZ["k"]
    M = Matrix(R, 1, 1, [R.gen()])
    gamma = QQ(gamma)
    a = gamma.numerator()
    b = gamma.denominator()
    k = lambda p, a=a, b=b: -(-a * p // b)
    m = lambda p, e=e: p**e

    ans = remainder_forest(M, m, k, kbase=1, indices=prime_range(n))
    return {p: ans[p][0, 0] for p in prime_range(n)}


def batch_harmonic(n, e, gamma, j, proj=False):
    r"""
    Amortized computation of harmonic sums.

    INPUT::

     - `n`: bound on primes
     - `e`: order of the series expansion
     - `gamma`: a rational number
     - `j`: a positive integer

    OUTPUT::
       A dict whose value at a prime `p` is the truncated harmonic sum
          \sum_{k=1}^{ceil(gamma*p)-1} k^{-j} mod p^e.

       If `proj` is True, instead return pairs (x,y) representing x,y mod p^e.
    """
    # Represent summation as a 2x2 matrix.
    R = ZZ["k"]
    y = R.gen()
    M = Matrix(R, [[y**j, 0], [1, y**j]])

    gamma = QQ(gamma)
    a, b = gamma.as_integer_ratio()
    k = lambda p, a=a, b=b: -(-a * p // b)
    m = lambda p, e=e: p**e
    V = Matrix(ZZ, [[0, 1]])
    repeated_entries = [-1,-1,-1,0]

    ans = remainder_forest(M, m, k, kbase=1, indices=prime_range(n), V=V, repeated_entries=repeated_entries)
    if proj:
        return ans
    return {p: mat[0, 0] * mat[0, 1].inverse_mod(p**e) for p, mat in ans.items()}
