# Kdirectional

R-package of tools for anisotropy analysis of point location data, also known as point patterns, based on the pairwise differences. Focus on tools that work both in 2D and in 3D.



## How to install

As usual, use `devtools` package:

```
library(devtools)
install_github('antiphon/Kdirectional')
```

To also install the documentation (vignette), use

```
install_github('antiphon/Kdirectional', build_vignettes = TRUE)
```

The package depends on the packages `ellipsoid` and `sphere`, which you can install with

```
install_github('antiphon/sphere')
install_github('antiphon/ellipsoid')
```
