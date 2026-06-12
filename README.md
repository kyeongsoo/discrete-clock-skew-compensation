# discrete-clock-skew-compensation

This repository provides the C implementation of the algorithms
for clock skew compensation based on integer linear scaling rounded
to the nearest integer and the results of the investigation of their
space-time trade-off [^1].

## Prerequites

- [The GNU MPFR library](https://www.mpfr.org/)
- [The GNU GMP library](https://gmplib.org/)

## Platforms

The code has been tested on Linux and Windows with [WSL](https://learn.microsoft.com/en-us/windows/wsl/) and [MSYS2](https://www.msys2.org/).

## References

[^1]: Kyeong Soo Kim, "Space-time trade-off in integer linear scaling rounded to the nearest integer through multiplicative and additive decomposition," arXiv e-prints [arXiv:2605.21400v](https://arxiv.org/abs/2605.21400) [cs.DS], May 2026.
