# pyrforest

This package is a simple wrapper for the [rforest](pyrforest/lib/README.md) library code into SageMath.

## Install



```
sage -pip install --upgrade  git+https://github.com/edgarcosta/pyrforest.git@master
```

If you don't have permissions to install it system wide, please add the flag ``--user`` to install it just for you.

```
sage -pip install --user --upgrade  git+https://github.com/edgarcosta/pyrforest.git@master
```




## Functions provided
### Interface with rforest
#### remainder_forest(M, m, k, kbase=0, indices=None, V=None, ans=None, kappa=None, projective=False)
#### remainder_forest_generic_prime(M, d, e, k, indices=None, m=None, kbase=0, V=None, ans=None, kappa=None)

### Examples of usage

- `batch_factorial(n, e, gamma)`

  Return a dict whose value at a prime $p$ equals $(\lceil \gamma p \rceil-1)! \pmod p^{e}$.
 
- `batch_harmonic(n, e, gamma, j, proj=False)`

  Return a dict whose value at a prime $p$ is the truncated harmonic sum
  
  $$ \sum_{k=1}^{\lceil \gamma p \rceil-1} k^{-j} mod p^e. $$
  
  If `proj` is True, instead return a $1 \times 2$ matrix $[x, y]$ representing $x/y \pmod{p^e}$.ƒ

