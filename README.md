# pyrforest

This package is a simple wrapper for the [rforest](https://github.com/edgarcosta/rforest) library code into SageMath.

## Install

```
sage -pip install --no-build-isolation --upgrade git+https://github.com/edgarcosta/pyrforest.git
```

If you don't have permissions to install it system wide, please add the flag ``--user`` to install it just for you.

```
sage -pip install --user --no-build-isolation --upgrade git+https://github.com/edgarcosta/pyrforest.git
```

## Development

Clone with submodules:
```
git clone --recurse-submodules https://github.com/edgarcosta/pyrforest.git
```

If you already cloned without `--recurse-submodules`:
```
git submodule update --init
```

Then install in editable mode:
```
make install
```

Run tests:
```
make test
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

  If `proj` is True, instead return a $1 \times 2$ matrix $[x, y]$ representing $x/y \pmod{p^e}$.
