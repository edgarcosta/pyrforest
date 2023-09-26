# pyrforest

This package is a simple wrapper for the [rforest](pyrforest/lib/README.md) library code into SageMath.


In particular, we provide:
- `remainder_forest(M, m, k, kbase=0, indices=None, V=None, ans=None, kappa=None, projective=False)`
- `remainder_forest_generic_prime(M, d, e, k, indices=None, m=None, kbase=0, V=None, ans=None, kappa=None)`

and two simple use cases:
- `batch_factorial(n, e, gamma)`
- `batch_harmonic(n, e, gamma, j, proj=False)`
Return dict whose value at a prime `p` is the truncated harmonic sum
          $$ \sum_{k=1}^{\lceilgamma*p \rceil-1} k^{-j} mod p^e. $$

       If `proj` is True, instead return pairs (x,y) representing x,y mod p^e.

