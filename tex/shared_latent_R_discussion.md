# Shared latent R discussion

Date: 2026-06-17

## Context

We started from the current SuSiE-style variational model with uncertainty in
the LD matrix `R`. The likelihood was updated to the summary-statistic form

```tex
x \mid \cdots \sim N(mean, \bar R / N).
```

The current implementation in `code/new-model.R` uses one independent latent
inverse-Wishart matrix per effect:

```tex
R^{(\ell)} \sim IW_J(\nu_0 \bar R, \nu_0 + J + 1).
```

Conditional on effect `l` selecting coordinate `j`, the selected column is
represented as

```tex
R^{(\ell)}_j = s_{\ell j} u_j(z_{\ell j}).
```

## Empirical behavior observed

Simulation experiments showed:

- With `L = 1`, `J = 100`, `N = 10000`, and `nu0_true = 1000`, the model
  recovered both the true coordinate and `nu0` well.
- With `L = 2`, coordinate recovery by aggregate PIP can remain good, but the
  estimated `nu0` is much too small or unstable.
- Contribution plots showed that a single effect can absorb both true signals.
  The estimated contribution for effect 1 can contain a fake peak at the second
  signal location.
- In fixed-`nu0` runs, the two estimated effect contributions can develop
  artificial peaks of opposite signs that cancel in the total fitted mean.
- Starting `nu0` large or warm-starting the fixed-`nu0` fit at the EB-recovered
  coordinates did not fix this behavior in the current CAVI scheme.

## Interpretation

The main issue appears to be model flexibility, not just initialization.

The variational update for the latent column shape has the form

```tex
m_z = \mu_j + \alpha_j (y_{-j} - \mu_j y_j),
\qquad
\alpha_j = \frac{N T_{1j}}{N T_{2j} + \rho_j E[\omega_j]}.
```

When `N` is large, `alpha_j` can be close to 1 even for large absolute values
of `nu0`. Thus each effect-specific latent column can move strongly toward the
current residual. This lets one effect invent residual-shaped peaks and lets
different effects create incompatible columns whose artifacts cancel.

The current model therefore treats uncertainty in `R` partly as an extra
flexible regression basis. This can improve the variational objective while
producing poor decomposition and poor `nu0` recovery.

## Shared latent R idea

A more coherent model is to use one shared latent LD matrix:

```tex
R \sim IW_J(\nu_0 \bar R, \nu_0 + J + 1),
```

and

```tex
x \mid R, \beta, \gamma
\sim
N\left(
  \sum_{\ell=1}^L \beta_\ell R_{\gamma_\ell},
  \bar R / N
\right).
```

This prevents each effect from inventing its own incompatible LD column. The
symmetry constraint `R_{jk} = R_{kj}` couples the selected columns.

Exact inference under a shared inverse-Wishart `R` is hard because columns of
`R` are dependent.

## Preferred approximate direction

We discussed replacing the exact inverse-Wishart column model with a Gaussian
perturbation model:

```tex
R = \bar R + E,
```

where `E` is a dense symmetric Gaussian perturbation. We are not trying to
estimate `R` itself; `E` is nuisance uncertainty to integrate out for robust
fine-mapping.

A simple version is

```tex
E_{ab} = E_{ba}, \qquad E_{ab} \sim N(0, \tau^2_{ab})
```

independently over `a <= b`. Initially, use a scalar uncertainty parameter:

```tex
\tau^2 = c / \nu_0
```

or

```tex
\tau^2 = 1 / \nu_0.
```

Then

```tex
x =
\sum_{\ell=1}^L \beta_\ell \bar R_{\gamma_\ell}
+
\sum_{\ell=1}^L \beta_\ell E_{\gamma_\ell}
+
\epsilon,
\qquad
\epsilon \sim N(0, \bar R / N).
```

Since `E` is Gaussian, it can be integrated out conditional on `beta` and
`gamma`:

```tex
x \mid \beta, \gamma, \nu_0
\sim
N\left(
  \sum_{\ell=1}^L \beta_\ell \bar R_{\gamma_\ell},
  \bar R / N
  +
  \operatorname{Var}_E\left[
    \sum_{\ell=1}^L \beta_\ell E_{\gamma_\ell}
  \right]
\right).
```

This treats LD uncertainty as extra covariance in the observation model rather
than as effect-specific flexible columns. It should reduce fake canceling peaks.

## Open design questions

- Should `E` use one scalar variance parameter or entry-specific variances?
- Should entry-specific variances be calibrated from inverse-Wishart moments?
- Do we need to enforce positive definiteness of `R = \bar R + E`, or is a
  local Gaussian approximation around `\bar R` sufficient?
- How should the induced covariance

  ```tex
  \operatorname{Var}_E\left[
    \sum_\ell \beta_\ell E_{\gamma_\ell}
  \right]
  ```

  be approximated under the SuSiE variational distribution?
- Can we implement a first version by iteratively inflating the residual
  covariance using current variational moments and then running SuSiE-style
  updates with columns of `\bar R`?

## Current working hypothesis

The independent per-effect inverse-Wishart column model is too flexible for
multi-effect fine-mapping. A shared dense Gaussian perturbation `E` is a more
promising direction because it preserves global LD consistency while allowing
LD uncertainty to be integrated out.

