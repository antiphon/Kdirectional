# Tools for directional analysis of 2- and 3-D point patterns

## Description

`Kdirectional` is an R-package containing tools for anisotropy analysis of point location data, also known as point patterns, based on pairwise differences. Focus is on tools that work both in 2D and in 3D. 

The package was born as a collection of custom functions for manipulating and visualising 3D point patterns, such as the ice sheet bubble data analysed in our geometric anistropy estimation paper Rajala et al. (2016)[^ice]. There might still be some exotic functions functions hanging around, undocumented, from early 2010s.

In Rajala, Redenbach, Särkkä, Sormani (2018)[^rev] we collected "classical" point pattern statistics for anisotropy analysis, and in the follow up paper [^rev2] compared some of them in some simulation trials to see which "features" are most informative of anisotropy, and most powerful (statistically speaking) in testing for isotropy. `Kdirectional` includes all the summaries discussed in those papers. 

The package has absorbed most of the functionality of `ellipsoid` (https://www.github.com/antiphon/ellipsoid) in order to reduce dependencies. 

Then name comes from the first function, Ripley's K-function for directional analysis (see `Kdirectional::Kest_directional`).

The codebase starts to be almost decade old in places, some documentation is missing, and other is incomplete. Many of the functions are not based on proper research, and some are not properly tested. If you spot a point of improvement please contact the author or submit an issue at the github repository, https://www.github.com/antiphon/Kdirectional.

## How to install

To install from Github, use `devtools` or `remotes`. 

```r
# Requirements:
require(spatstat)
require(sp)
require(rgl)
require(mvtnorm)

# Install:
library(devtools) # or remotes
install_github('antiphon/Kdirectional')
```

To also install the documentation (vignette), use

```r
install_github('antiphon/Kdirectional', build_vignettes = TRUE)
```

The documentation should be available also as a `pkgdown` site.


[^ice]: T. Rajala, A. Särkkä, C. Redenbach, and M. Sormani. ‘Estimating Geometric Anisotropy in Spatial Point Patterns’. Spatial Statistics 15 (February 2016): 100–114. https://doi.org/10.1016/j.spasta.2015.12.005.

[^rev]: T. Rajala, C. Redenbach, A. Särkkä, and M. Sormani. ‘A Review on Anisotropy Analysis of Spatial Point Patterns’. Spatial Statistics 28 (December 2018): 141–68. https://doi.org/10.1016/j.spasta.2018.04.005.

[^rev2]: T. Rajala, C. Redenbach, A. Särkkä, and M. Sormani. ‘Tests for Isotropy in Spatial Point Patterns – A Comparison of Statistical Indices’. Spatial Statistics 52 (December 2022): 100716. https://doi.org/10.1016/j.spasta.2022.100716.
