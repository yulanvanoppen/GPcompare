# GPcompare
_Compare groups of time series using Guassian processes (GPs) through posterior mean difference distribution evaluation._

&nbsp;


## Contents
`[p, T, E] = GPcompare(X, Y, ...)` compares two groups of time series data given in 2-D arrays `X` and `Y` with rows corresponding to (aligned) time points and columns to individuals (or repetitions). The dimensions of `X` and `Y` are not restricted to be equal.

**Output:**
- `p` tail probability of the zero vector;
- `T` effect size;
- `E` effect size decomposition such that T = E'E.

&nbsp;


## Example
Run `example.m` to reproduce the results given in Figure 1 and 2 below.

<img src="https://github.com/yulanvanoppen/GPcompare/blob/master/figures/minimal_difference.svg" width="800">

**Figure 1.** Synthetic time series data from two populations generated from non-stationary Gaussian processes (left and middle panels). The two populations, **X** and **Y**, are assumed to differ slightly in the parameterization of the deterministic part. Sample single-cell trajectories generated from each population are shown as thin, solid lines, along with the inferred population (posterior) mean distributions, which are summarized by their means (thick, solid) ± 2 standard deviations (dashed). The right panel depicts the distribution of the difference of population (posterior) means, summarized by its mean (solid) ± 2 standard deviations (dashed), while the bars represent the decomposed effect size across a sparse grid of time points. Note that the bands around the solid line include zero for most time points, and the effect size components are relatively small. The tail probability for this comparison is 𝜀 = 0.116 with an effect size of T = 14.2, confirming that the two population means do not seem to differ too much.

<img src="https://github.com/yulanvanoppen/GPcompare/blob/master/figures/observable_difference.svg" width="800">

**Figure 2.** A similar setup as in Figure 1, but with the **Y** data generated using a significantly perturbed mean in the second half of the cell cycle. This time, the bands around the solid line on the right panel do not include zero for a large fraction of the time points, and the effect size components are substantially larger in the second half of the cell cycle compared to the first (and also compared to the effect sizes of the corresponding panel of the first row). The tail probability for this comparison is 𝜀 = 2.92E-16 with an effect size of T = 96.4, suggesting the population means are unlikely to be equal.

&nbsp;


## Technical details

See Methods in [1].

&nbsp;


## References

[1] Guerra, Paolo, Vuillemenot, Luc-Alban P.E., van Oppen, Yulan B., Been, Marije, & Milias-Argeitis, Andreas (in press). TORC1 and PKA activity towards ribosome biogenesis oscillates in synchrony with the budding yeast cell cycle. _Journal of Cell Science_.

&nbsp;


## DISCLAIMER
_All software in this repository is covered by this disclaimer:_

While every effort is made to deliver high quality products, no guarantee is made that the products are free from defects. The software is provided "as is", and you use the software at your own risk.

No warranties are made as to performance, merchantability, fitness for a particular purpose, or any
other warranties whether expressed or implied.

No oral or written communication from or information provided by the author shall create a warranty.

Under no circumstances shall the author be liable for direct, indirect, special, incidental, or consequential damages resulting from the use, misuse, or inability to use this software, even if the author has been advised of the possibility of such damages.
