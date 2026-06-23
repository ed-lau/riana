# Robust CV for turnover rate constants — basis for `1.4826·MAD(ln k)`

- **Date:** 2026-06-23
- **Status:** reference — derivation + primary citations; estimator sound, **bias < 2.5% for σ_log ≲ 0.30**
- **Scope:** per-curve / per-precursor dispersion of the fitted first-order turnover rate constant *k* (turnover fitting-error reporting)
- **Code:** turnover fit (`k_deg`) in [`riana/core/fitting.py`](../riana/core/fitting.py); the within-protein geometric robust CV `1.4826·MAD(ln k)` is computed in [`tests/benchmark/bench_mbr_ab.py`](../tests/benchmark/bench_mbr_ab.py) (`kcv` inside `_protein_metrics`) — per `(protein id, condition)` over peptide-level `k>0`, median over proteins. **Verified 2026-06-23: that formula matches this derivation exactly** (constant, MAD-of-log, k>0 filter).
- **Verify:** [§Reproduce](#reproduce) regenerates the constant, the bias table, and the robustness demo from scratch

## Question

We summarise the dispersion of a fitted turnover rate constant *k* as a **geometric
robust CV**:

$$\widehat{\mathrm{CV}}_{\text{geo,robust}} = 1.4826 \cdot \mathrm{MAD}(\ln k)$$

over a set of *k* estimates (replicate fits, bootstrap resamples, or peptide-level
estimates within a protein). Three things need to be on the record before this is
relied on downstream:

1. **What is it actually estimating?** Is `1.4826·MAD(ln k)` a CV in any rigorous
   sense, or just a convenient proxy?
2. **Is there a single citation, or is it a composition?**
3. **Over what range is it valid**, and where do we switch to something exact?

## TL;DR

- **It is a composition of two standard results, not one named formula.** No single
  paper writes `1.4826·MAD(ln k)`; it is *(robust scale on the log scale)* combined
  with *(the log-scale SD ≈ the geometric CV)*. Each half has solid primary
  citations.
- **`1.4826·MAD(ln k)` is a Fisher-consistent, outlier-resistant estimate of σ_log**,
  the standard deviation of `ln k`. The multiplier is exactly
  $1/\Phi^{-1}(3/4)=1.4826\ldots$ (Hampel 1974; Rousseeuw & Croux 1993).
- **σ_log is the small-dispersion limit of the geometric CV** of a lognormal *k*: the
  exact CV is $\sqrt{e^{\sigma^2}-1}$, and $\sqrt{e^{\sigma^2}-1}\to\sigma$ as
  σ→0 (Limpert et al. 2001).
- **So reporting it as a CV is exact to first order.** Verified bias (it always
  *under*-reports): **−0.25% at σ_log = 0.10, −1.0% at 0.20, −2.2% at 0.30, −6.2% at
  0.50.**
- **Robust is the point.** MAD has a 50% breakdown point — a handful of badly-fit
  peptides cannot blow up the dispersion the way they do for SD/mean. In a quick
  check, one gross outlier in 400 lognormal draws (true σ = 0.20) pushes naive
  `SD(ln k)` to **0.40** while the robust estimate holds at **0.20**.
- **Recommendation:** keep `1.4826·MAD(ln k)` for well-behaved curves
  (σ_log ≲ 0.30, i.e. CV ≲ 30%); past that, report the exact
  $\sqrt{e^{\hat\sigma^2}-1}$ with $\hat\sigma = 1.4826\cdot\mathrm{MAD}(\ln k)$.
  One-line change — see [§Reproduce](#reproduce).

## Setup & notation

Given a set of rate-constant estimates $\{k_1,\dots,k_n\}$ with $k_i>0$, let
$y_i=\ln k_i$.

- **Lognormal assumption.** We treat *k* as lognormal, i.e. $\ln k \sim
  \mathcal N(\mu,\sigma^2)$. Justification for turnover constants: *k* is strictly
  positive, fit/measurement error on rates is multiplicative rather than additive,
  and fold-changes are symmetric in log space. This is the same assumption that
  motivates reporting *geometric* means for *k*.
- **MAD** $= \operatorname{med}_i\,\lvert y_i - \operatorname{med}_j y_j\rvert$ — the
  median of absolute deviations from the median.
- **σ_log** — the SD of $\ln k$; the dispersion parameter we estimate.
- **GCV** (geometric CV) — the coefficient of variation of the lognormal *k* on its
  natural scale $= \sqrt{e^{\sigma^2}-1}$.

## Derivation

### 1. `1.4826·MAD` is a robust estimate of σ_log

For $Y\sim\mathcal N(\mu,\sigma^2)$ the population MAD satisfies
$P(\lvert Y-\mu\rvert \le \mathrm{MAD})=\tfrac12$. Standardising,

$$\frac{\mathrm{MAD}}{\sigma}=\Phi^{-1}(0.75)\approx 0.674490
\quad\Longrightarrow\quad
\hat\sigma_{\log}=\frac{\mathrm{MAD}(\ln k)}{\Phi^{-1}(0.75)} = 1.4826\cdot\mathrm{MAD}(\ln k).$$

The multiplier is exactly $1/\Phi^{-1}(3/4)=1/0.674490\approx 1.482602$ — the
reciprocal of the standard-normal upper quartile — chosen so the scaled MAD is
**Fisher-consistent** for σ at the normal model. MAD as a robust measure of scale
traces to Hampel (1974); its status as the reference robust scale estimator (50%
breakdown point, bounded influence) is established in Rousseeuw & Croux (1993). For
small *n* the asymptotic 1.4826 should carry a finite-sample correction $c_n$
(Croux & Rousseeuw 1992) — worth applying if curves are summarised from few
replicates.

### 2. The log-scale SD is the geometric CV in the small-dispersion limit

For lognormal *k* with $\ln k\sim\mathcal N(\mu,\sigma^2)$, the natural-scale moments
give a coefficient of variation that is independent of μ:

$$\mathrm{CV}[k]=\frac{\mathrm{SD}[k]}{\mathbb E[k]}=\sqrt{e^{\sigma^2}-1}\;\equiv\;\mathrm{GCV}.$$

Taylor-expanding for small σ,

$$\sqrt{e^{\sigma^2}-1}=\sigma\sqrt{1+\tfrac{\sigma^2}{2}+\cdots}\;\approx\;\sigma.$$

Hence $1.4826\cdot\mathrm{MAD}(\ln k)\approx\hat\sigma_{\log}\approx\mathrm{GCV}$.
The estimator *is* a CV to first order; the gap is the higher-order lognormal
curvature quantified in [§Validity bounds](#validity-bounds). The accepted
definition $\mathrm{GCV}=\sqrt{e^{\sigma^2}-1}$ — as opposed to Kirkwood's
deprecated $e^{\sigma}-1$ — is standard in the lognormal and pharmacometrics
literature (Limpert et al. 2001; it is also what SAS and NONMEM report).

### 3. Why it deserves the name "CV" (intuition)

Because the median commutes with the (monotonic) log,

$$\mathrm{MAD}(\ln k)=\operatorname{med}_i\bigl\lvert \ln k_i-\ln k_{\text{med}}\bigr\rvert
=\operatorname{med}_i\Bigl\lvert \ln\tfrac{k_i}{k_{\text{med}}}\Bigr\rvert
\approx\operatorname{med}_i\Bigl\lvert \tfrac{k_i-k_{\text{med}}}{k_{\text{med}}}\Bigr\rvert,$$

i.e. to first order it is a **robust median relative deviation** of *k* about its
median — exactly the quantity a CV is meant to capture, with the median playing the
role of the mean and the MAD the role of the SD. This is the log-scale analogue of
the raw-scale robust CV $\mathrm{RCV}_M = 1.4826\cdot\mathrm{MAD}/\operatorname{med}(k)$
used in chemometrics (Reimann et al. 2008; Varmuza & Filzmoser 2009; analysed in
Arachchige et al. 2022); the two agree to first order via the identity above.

## Validity bounds

`1.4826·MAD(ln k)` always **under-reports** the true geometric CV, because
$\sigma < \sqrt{e^{\sigma^2}-1}$. The bias is set by σ_log alone — independent of the
rate, of *n*, and of the median:

| σ_log = `1.4826·MAD(ln k)` | exact GCV = √(e^σ²−1) | rel. error of reporting σ_log as CV |
|---|---|---|
| 0.10 (≈10% CV) | 0.1003 | −0.25% |
| 0.20 (≈20% CV) | 0.2020 | −1.0% |
| 0.30 (≈30% CV) | 0.3069 | −2.2% |
| 0.50 (≈50% CV) | 0.5329 | −6.2% |

(Generated in [§Reproduce](#reproduce).) **Practical rule:** for well-fit turnover
curves with σ_log ≲ 0.30 the approximation is good to ~2% — below the noise floor of
the fits — so reporting `1.4826·MAD(ln k)` directly as the CV is fine. For noisier
peptides past that, substitute $\hat\sigma=1.4826\cdot\mathrm{MAD}(\ln k)$ into
$\sqrt{e^{\hat\sigma^2}-1}$ and report the exact GCV. The robust estimate of σ_log is
unchanged either way; only the final transform differs.

**Implementation status (2026-06-23).** `bench_mbr_ab.py` reports the first-order
`σ_log = 1.4826·MAD(ln k)` directly as the CV at *all* magnitudes — it does **not**
yet apply the exact `√(e^σ²−1)` transform past σ_log ≈ 0.30. This is fine for the
numbers we rely on (within-protein curated ~17%, R²>0.8 ~25% — both inside the
≲2.5%-bias band), but the **all-converged** population (~46–49% on LVE, the noisy
tail) is past the threshold and is therefore *under*-reported by ~5%; read those
large values as a robust spread descriptor, and apply the one-line exact transform if
a precise CV is needed there.

**One caveat on the assumption, not the algebra.** The 1.4826 factor and the
GCV ≈ σ_log link both assume approximate log-normality. MAD stays a sensible scale
measure under mild contamination and slight skew, but if the per-curve *k*
distribution is strongly multimodal or heavy-tailed beyond a few outliers, neither
"σ_log" nor "GCV" has a clean interpretation — treat the number as a robust spread
descriptor, not a distributional CV.

## Citations

**Median absolute deviation & the 1.4826 scale factor**

- Hampel, F. R. (1974). "The influence curve and its role in robust estimation."
  *Journal of the American Statistical Association*, 69(346), 383–393. — MAD as a
  robust measure of scale.
- Rousseeuw, P. J., & Croux, C. (1993). "Alternatives to the median absolute
  deviation." *Journal of the American Statistical Association*, 88(424),
  1273–1283. — MAD as the reference robust scale estimator; 50% breakdown point; the
  $1/\Phi^{-1}(3/4)$ consistency factor.
- Croux, C., & Rousseeuw, P. J. (1992). "Time-efficient algorithms for two highly
  robust estimators of scale." In Y. Dodge & J. Whittaker (Eds.), *Computational
  Statistics*, Vol. 1, 411–428. Physica-Verlag. — finite-sample correction factors
  $c_n$ for the MAD scale estimate.
- Huber, P. J. (1981). *Robust Statistics*. Wiley. — standard textbook treatment of
  MAD-based scale estimation.

**Robust CV via MAD (raw-scale analogue)**

- Arachchige, C. N. P. G., Prendergast, L. A., & Staudte, R. G. (2022). "Robust
  analogs to the coefficient of variation." *Journal of Applied Statistics*, 49(2),
  268–290. DOI: 10.1080/02664763.2020.1808599 (arXiv:1907.01110). — MAD/median and
  IQR/median as robust CVs; influence functions, bias, interval estimators.
- Reimann, C., Filzmoser, P., Garrett, R. G., & Dutter, R. (2008). *Statistical Data
  Analysis Explained: Applied Environmental Statistics with R*. Wiley. — robust CV
  via MAD in applied practice.
- Varmuza, K., & Filzmoser, P. (2009). *Introduction to Multivariate Statistical
  Analysis in Chemometrics*. CRC Press. — same, chemometrics context.

**Lognormal distribution & the geometric CV**

- Limpert, E., Stahel, W. A., & Abbt, M. (2001). "Log-normal distributions across
  the sciences: keys and clues." *BioScience*, 51(5), 341–352. — geometric
  mean / SD / CV; the $\sqrt{e^{\sigma^2}-1}$ definition.
- Aitchison, J., & Brown, J. A. C. (1957). *The Lognormal Distribution*. Cambridge
  University Press. — primary derivation of the lognormal moments and CV.
- Koopmans, L. H., Owen, D. B., & Rosenblatt, J. I. (1964). "Confidence intervals for
  the coefficient of variation for the normal and log normal distributions."
  *Biometrika*, 51(1–2), 25–32. DOI: 10.1093/biomet/51.1-2.25 — CV of the lognormal
  and its confidence intervals.
- Kirkwood, T. B. L. (1979). "Geometric means and measures of dispersion."
  *Biometrics*, 35(4), 908–909. — proposes the alternative GCV $= e^{\sigma}-1$;
  cited for completeness (this is *not* the definition used here).

## Reproduce

```python
# regenerates the consistency constant, the §Validity bias table, and the
# robustness demo. requires numpy + scipy.
import numpy as np
from scipy.stats import norm

c = 1.0 / norm.ppf(0.75)          # Fisher-consistency factor
assert abs(c - 1.4826022185) < 1e-9
print(f"1/Phi^-1(0.75) = {c:.10f}")

print("sigma_log | exact GCV sqrt(e^s^2-1) | rel.err")
for s in (0.10, 0.20, 0.30, 0.50):
    gcv = np.sqrt(np.exp(s**2) - 1.0)
    print(f"  {s:.2f}   |   {gcv:.4f}   |   {(s - gcv) / gcv * 100:+.2f}%")

def geo_robust_cv(k, exact=False):
    """Robust geometric CV of positive rate constants k.
    exact=False -> 1.4826*MAD(ln k)         (sigma_log; first-order CV)
    exact=True  -> sqrt(exp(sigma^2) - 1)   (exact lognormal CV)
    """
    y = np.log(np.asarray(k, float))
    sigma = 1.482602218505602 * np.median(np.abs(y - np.median(y)))
    return np.sqrt(np.exp(sigma ** 2) - 1.0) if exact else sigma

# robustness demo: 400 lognormal draws at true sigma=0.20, plus one gross outlier
rng = np.random.default_rng(0)
k = np.exp(rng.normal(np.log(0.05), 0.2, size=400))
k[0] = 50.0
print()
print(f"robust  sigma_log estimate : {geo_robust_cv(k):.4f}  (true 0.20)")
print(f"robust  exact GCV          : {geo_robust_cv(k, exact=True):.4f}")
print(f"naive   SD(ln k)           : {np.std(np.log(k)):.4f}  <- inflated by the outlier")
```

Expected output:

```
1/Phi^-1(0.75) = 1.4826022185
sigma_log | exact GCV sqrt(e^s^2-1) | rel.err
  0.10   |   0.1003   |   -0.25%
  0.20   |   0.2020   |   -1.00%
  0.30   |   0.3069   |   -2.24%
  0.50   |   0.5329   |   -6.18%

robust  sigma_log estimate : 0.2035  (true 0.20)
robust  exact GCV          : 0.2056
naive   SD(ln k)           : 0.3986  <- inflated by the outlier
```