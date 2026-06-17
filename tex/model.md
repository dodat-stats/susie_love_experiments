::: center
**susie withh uncertainty in R**
:::

# Single-effect regression with uncertainty in $R$

Let $x \in \mathbb R^J$ denote the observed vector of summary
statistics, computed from a study with sample size $N$. We consider the
single-effect model
$$x \mid R,\beta,\gamma=j
\sim
N_J\left(\beta R_j,\frac{\bar R}{N}\right),$$ where $R_j$ is the $j$-th
column of $R$,
$$\beta \sim N(0,\sigma_0^2),
\qquad
\gamma \sim \operatorname{Unif}\{1,\ldots,J\},$$ and
$$R \sim IW_J(\nu_0 \bar R,\nu_0+J+1).$$ Here $\bar R$ is treated as
known, while $\nu_0$ controls the amount of uncertainty around $\bar R$.
Larger $\nu_0$ corresponds to stronger concentration of $R$ around
$\bar R$.

Let $$K=\bar R^{-1}.$$ The likelihood precision is therefore $N K$. For
each candidate index $j$, define
$$\rho_j=\bar R_{jj},$$ $$\mu_j=\frac{\bar R_{-j,j}}{\bar R_{jj}},$$ and
$$Q_j
=
\bar R_{-j,-j}
-
\frac{\bar R_{-j,j}\bar R_{j,-j}}{\bar R_{jj}}.$$

For the $j$-th column of $R$, write $$s=R_{jj},
\qquad
z=\frac{R_{-j,j}}{R_{jj}}.$$ Then $$R_j=s u_j(z),$$ where
$u_j(z)\in \mathbb R^J$ is the vector whose $j$-th entry is $1$ and
whose remaining entries are $z$. That is, $$u_j(z)_j=1,
\qquad
u_j(z)_{-j}=z.$$

Under the inverse-Wishart prior, $$s\mid \nu_0,\gamma=j
\sim
\operatorname{InvGamma}(a,b_j),$$ where $$a=\frac{\nu_0+2}{2},
\qquad
b_j=\frac{\nu_0\rho_j}{2}.$$ Also, $$z\mid \nu_0,\gamma=j
\sim
t_{\nu_0+3}
\left(
\mu_j,
\frac{Q_j}{(\nu_0+3)\rho_j}
\right).$$

We use the exact scale-mixture representation of the multivariate
$t$-distribution. Introduce $$\omega\mid \nu_0
\sim
\chi^2_{\nu_0+3}.$$ Equivalently, $$\omega\mid \nu_0
\sim
\operatorname{Gamma}
\left(
c,\frac12
\right),
\qquad
c=\frac{\nu_0+3}{2},$$ using the shape-rate parameterization.
Conditional on $\omega$, $$z\mid \omega,\gamma=j
\sim
N_{J-1}
\left(
\mu_j,
\frac{Q_j}{\omega\rho_j}
\right).$$

Thus, for fixed $\gamma=j$, the augmented model is
$$x\mid \beta,s,z,\gamma=j
\sim
N_J\left(\beta s u_j(z),\frac{\bar R}{N}\right),$$
$$\beta\sim N(0,\sigma_0^2),$$
$$s\mid \nu_0,j\sim \operatorname{InvGamma}(a,b_j),$$
$$\omega\mid \nu_0\sim \operatorname{Gamma}\left(c,\frac12\right),$$ and
$$z\mid \omega,j\sim
N_{J-1}
\left(
\mu_j,
\frac{Q_j}{\omega\rho_j}
\right).$$

## Variational family

We use the variational approximation $$q(\gamma,\beta,s,\omega,z)
=
q(\gamma)q(\beta,s,\omega,z\mid \gamma).$$ Let $$q(\gamma=j)=\pi_j.$$
Conditional on $\gamma=j$, use the mean-field factorization
$$q_j(\beta,s,\omega,z)
=
q_j(\beta)q_j(s)q_j(\omega)q_j(z).$$ We take
$$q_j(\beta)=N(m_{\beta j},v_{\beta j}),$$
$$q_j(\omega)=\operatorname{Gamma}(\alpha_{\omega j},\beta_{\omega j}),$$
and $$q_j(z)=N_{J-1}(m_{zj},V_{zj}).$$ The factor $q_j(s)$ is a
one-dimensional nonstandard distribution on $(0,\infty)$, handled by
quadrature.

Define the variational moments $$\bar\beta_j=E_j[\beta]=m_{\beta j},$$
$$\overline{\beta^2}_j=E_j[\beta^2]=v_{\beta j}+m_{\beta j}^2,$$
$$\bar s_j=E_j[s],
\qquad
\overline{s^2}_j=E_j[s^2],$$ and $$\bar\omega_j=E_j[\omega].$$

For a working response $y$, define $$d_j(y)
=
y^TKy-\frac{y_j^2}{\rho_j}.$$ Equivalently, $$d_j(y)
=
(y_{-j}-\mu_j y_j)^TQ_j^{-1}(y_{-j}-\mu_j y_j).$$

Define $$T_{1j}
=
\bar\beta_j\bar s_j,$$ and $$T_{2j}
=
\overline{\beta^2}_j\overline{s^2}_j.$$ Then set $$C_{zj}
=
N T_{2j}+\rho_j\bar\omega_j.$$

## Coordinate ascent updates

The optimal Gaussian update for $q_j(z)$ is $$V_{zj}
=
\frac{Q_j}{C_{zj}},$$ and $$m_{zj}
=
\mu_j
+
\frac{N T_{1j}}{C_{zj}}
\left(
y_{-j}-\mu_j y_j
\right).$$ Define $$\alpha_j=\frac{N T_{1j}}{C_{zj}}.$$

Two scalar summaries are useful: $$B_j
=
E_j\left[y^T K u_j(z)\right],$$ and $$A_j
=
E_j\left[u_j(z)^T K u_j(z)\right].$$ Using the above Gaussian update,
$$B_j
=
\frac{y_j}{\rho_j}
+
\alpha_j d_j(y),$$ and $$A_j
=
\frac1{\rho_j}
+
\alpha_j^2 d_j(y)
+
\frac{J-1}{C_{zj}}.$$

The optimal Gaussian update for $q_j(\beta)$ is $$v_{\beta j}^{-1}
=
\frac1{\sigma_0^2}
+
N\overline{s^2}_j A_j,$$ and $$m_{\beta j}
=
v_{\beta j}N\bar s_j B_j.$$

The optimal one-dimensional update for $q_j(s)$ is $$q_j(s)
\propto
s^{-(a+1)}
\exp
\left\{
-\frac{b_j}{s}
+
D_j s
-
\frac12 H_j s^2
\right\},
\qquad s>0,$$ where $$D_j=N m_{\beta j}B_j,$$ and $$H_j=
N\overline{\beta^2}_j A_j.$$ Equivalently, using the transformation
$s=e^\xi$, $$\log q_j(\xi)
=
\text{constant}
-
a\xi
-
b_j e^{-\xi}
+
D_j e^\xi
-
\frac12 H_j e^{2\xi}.$$ The moments $$E_j[s],
\qquad
E_j[s^2],
\qquad
E_j[\log s],
\qquad
E_j[1/s]$$ can be computed by one-dimensional quadrature over $\xi$.

Next define $$\Delta_j
=
\rho_j
E_j\left[
(z-\mu_j)^TQ_j^{-1}(z-\mu_j)
\right].$$ Using the Gaussian update for $q_j(z)$, $$\Delta_j
=
\rho_j
\left[
\alpha_j^2 d_j(y)
+
\frac{J-1}{C_{zj}}
\right].$$ Then the optimal update for $q_j(\omega)$ is $$q_j(\omega)
=
\operatorname{Gamma}(\alpha_{\omega j},\beta_{\omega j}),$$ with
$$\alpha_{\omega j}
=
\frac{\nu_0+J+2}{2},$$ and $$\beta_{\omega j}
=
\frac{1+\Delta_j}{2}.$$ Therefore, $$E_j[\omega]
=
\frac{\alpha_{\omega j}}{\beta_{\omega j}},$$ and $$E_j[\log\omega]
=
\psi(\alpha_{\omega j})-\log\beta_{\omega j}.$$

The assignment probabilities are updated by $$\pi_j
=
\frac{\exp\{\mathcal F_j\}}
{\sum_{k=1}^J \exp\{\mathcal F_k\}},$$ where $\mathcal F_j$ is the local
variational lower bound contribution for candidate $j$: $$\mathcal F_j
=
E_j
\left[
\log p_{\nu_0}(s\mid j)
+
\log p_{\nu_0}(\omega)
+
\log p(z\mid \omega,j)
+
\log p(\beta)
+
\log p(y\mid \beta,s,z,j)
-
\log q_j(s)
-
\log q_j(\omega)
-
\log q_j(z)
-
\log q_j(\beta)
\right].$$

## Variational empirical Bayes update for $\nu_0$

For fixed variational distribution $q$, update $\nu_0$ by maximizing the
part of the ELBO depending on $\nu_0$. Since the $z\mid\omega,j$ density
does not depend on $\nu_0$ under the latent-scale representation, the
relevant objective is $$\mathcal M(\nu_0)
=
\sum_{j=1}^J
\pi_j
\left[
E_j\log p_{\nu_0}(s\mid j)
+
E_j\log p_{\nu_0}(\omega)
\right].$$ Here $$s\mid \nu_0,j\sim \operatorname{InvGamma}(a,b_j),
\qquad
a=\frac{\nu_0+2}{2},
\qquad
b_j=\frac{\nu_0\rho_j}{2},$$ so $$E_j\log p_{\nu_0}(s\mid j)
=
a\log b_j
-
\log\Gamma(a)
-
(a+1)E_j[\log s]
-
b_jE_j[1/s].$$ Also $$\omega\mid \nu_0
\sim
\operatorname{Gamma}\left(c,\frac12\right),
\qquad
c=\frac{\nu_0+3}{2},$$ so $$E_j\log p_{\nu_0}(\omega)
=
c\log\frac12
-
\log\Gamma(c)
+
(c-1)E_j[\log\omega]
-
\frac12 E_j[\omega].$$ Thus the empirical Bayes update is
$$\nu_0^{\text{new}}
=
\arg\max_{\nu_0}
\mathcal M(\nu_0).$$ If finite inverse-Wishart variance is desired,
impose $\nu_0>2$ by parameterizing $$\nu_0=2+\exp(\eta).$$

# SuSiE with uncertainty in $R$

We now extend the single-effect model to a sum of $L$ single effects.
For each effect $\ell=1,\ldots,L$,
$$\gamma_\ell\sim \operatorname{Unif}\{1,\ldots,J\},$$
$$\beta_\ell\sim N(0,\sigma_\ell^2),$$ and independently $$R^{(\ell)}
\sim
IW_J(\nu_0\bar R,\nu_0+J+1).$$ Conditional on $\gamma_\ell=j$, the
$\ell$-th effect uses the $j$-th column $$R^{(\ell)}_j.$$ The
observation model is $$x
\mid
\{\beta_\ell,\gamma_\ell,R^{(\ell)}\}_{\ell=1}^L
\sim
N_J
\left(
\sum_{\ell=1}^L
\beta_\ell R^{(\ell)}_{\gamma_\ell},
\frac{\bar R}{N}
\right).$$

For each effect $\ell$ and candidate index $j$, define $$s_{\ell j}
=
R^{(\ell)}_{jj},
\qquad
z_{\ell j}
=
\frac{R^{(\ell)}_{-j,j}}{R^{(\ell)}_{jj}}.$$ Then $$R^{(\ell)}_j
=
s_{\ell j}u_j(z_{\ell j}).$$ As in the single-effect case,
$$s_{\ell j}\mid \nu_0,\gamma_\ell=j
\sim
\operatorname{InvGamma}(a,b_j),$$ with $$a=\frac{\nu_0+2}{2},
\qquad
b_j=\frac{\nu_0\rho_j}{2},$$ and $$\omega_{\ell j}\mid \nu_0
\sim
\operatorname{Gamma}
\left(
c,\frac12
\right),
\qquad
c=\frac{\nu_0+3}{2},$$ with
$$z_{\ell j}\mid \omega_{\ell j},\gamma_\ell=j
\sim
N_{J-1}
\left(
\mu_j,
\frac{Q_j}{\omega_{\ell j}\rho_j}
\right).$$

## Variational family

Use the SuSiE-style factorization $$q
=
\prod_{\ell=1}^L
q_\ell(\gamma_\ell,\beta_\ell,s_\ell,\omega_\ell,z_\ell).$$ For each
effect, $$q_\ell(\gamma_\ell=j)=\pi_{\ell j}.$$ Conditional on
$\gamma_\ell=j$, use
$$q_{\ell j}(\beta_\ell,s_{\ell j},\omega_{\ell j},z_{\ell j})
=
q_{\ell j}(\beta_\ell)
q_{\ell j}(s_{\ell j})
q_{\ell j}(\omega_{\ell j})
q_{\ell j}(z_{\ell j}),$$ where $$q_{\ell j}(\beta_\ell)
=
N(m_{\beta,\ell j},v_{\beta,\ell j}),$$ $$q_{\ell j}(\omega_{\ell j})
=
\operatorname{Gamma}(\alpha_{\omega,\ell j},\beta_{\omega,\ell j}),$$
and $$q_{\ell j}(z_{\ell j})
=
N_{J-1}(m_{z,\ell j},V_{z,\ell j}).$$ The factor
$q_{\ell j}(s_{\ell j})$ is one-dimensional and updated by quadrature.

## SuSiE-style coordinate ascent updates

Let $$m_\ell
=
E_q[\beta_\ell R^{(\ell)}_{\gamma_\ell}]$$ be the current posterior mean
contribution from effect $\ell$. Then $$m_\ell
=
\sum_{j=1}^J
\pi_{\ell j}
E_{\ell j}[\beta_\ell]
E_{\ell j}[s_{\ell j}]
E_{\ell j}[u_j(z_{\ell j})].$$

To update effect $\ell$, define the partial residual $$y_\ell
=
x-\sum_{k\neq \ell}m_k.$$ Then the update for effect $\ell$ is identical
to the single-effect update with response $y_\ell$ and prior variance
$\sigma_\ell^2$.

For each candidate $j$, define $$d_j(y_\ell)
=
y_\ell^T K y_\ell
-
\frac{y_{\ell j}^2}{\rho_j}.$$ Let $$\bar\beta_{\ell j}
=
E_{\ell j}[\beta_\ell],
\qquad
\overline{\beta^2}_{\ell j}
=
E_{\ell j}[\beta_\ell^2],$$ $$\bar s_{\ell j}
=
E_{\ell j}[s_{\ell j}],
\qquad
\overline{s^2}_{\ell j}
=
E_{\ell j}[s_{\ell j}^2],$$ and $$\bar\omega_{\ell j}
=
E_{\ell j}[\omega_{\ell j}].$$ Set $$T_{1,\ell j}
=
\bar\beta_{\ell j}\bar s_{\ell j},$$ $$T_{2,\ell j}
=
\overline{\beta^2}_{\ell j}\overline{s^2}_{\ell j},$$ and $$C_{z,\ell j}
=
N T_{2,\ell j}
+
\rho_j\bar\omega_{\ell j}.$$ Then $$V_{z,\ell j}
=
\frac{Q_j}{C_{z,\ell j}},$$ and $$m_{z,\ell j}
=
\mu_j
+
\frac{N T_{1,\ell j}}{C_{z,\ell j}}
\left(
y_{\ell,-j}-\mu_j y_{\ell j}
\right).$$ Define $$\alpha_{\ell j}
=
\frac{N T_{1,\ell j}}{C_{z,\ell j}}.$$ Then $$B_{\ell j}
=
\frac{y_{\ell j}}{\rho_j}
+
\alpha_{\ell j}d_j(y_\ell),$$ and $$A_{\ell j}
=
\frac1{\rho_j}
+
\alpha_{\ell j}^2d_j(y_\ell)
+
\frac{J-1}{C_{z,\ell j}}.$$

The Gaussian update for $q_{\ell j}(\beta_\ell)$ is
$$v_{\beta,\ell j}^{-1}
=
\frac1{\sigma_\ell^2}
+
N\overline{s^2}_{\ell j}A_{\ell j},$$ and $$m_{\beta,\ell j}
=
v_{\beta,\ell j}
N
\bar s_{\ell j}
B_{\ell j}.$$

The one-dimensional update for $q_{\ell j}(s_{\ell j})$ is
$$q_{\ell j}(s)
\propto
s^{-(a+1)}
\exp
\left\{
-\frac{b_j}{s}
+
D_{\ell j}s
-
\frac12H_{\ell j}s^2
\right\},
\qquad s>0,$$ where $$D_{\ell j}
=
N m_{\beta,\ell j}B_{\ell j},$$ and $$H_{\ell j}
=
N
\left(
v_{\beta,\ell j}+m_{\beta,\ell j}^2
\right)
A_{\ell j}.$$ Using $s=e^\xi$, $$\log q_{\ell j}(\xi)
=
\text{constant}
-
a\xi
-
b_j e^{-\xi}
+
D_{\ell j}e^\xi
-
\frac12H_{\ell j}e^{2\xi}.$$ Moments of $s_{\ell j}$ are computed by
one-dimensional quadrature.

Next define $$\Delta_{\ell j}
=
\rho_j
E_{\ell j}
\left[
(z_{\ell j}-\mu_j)^TQ_j^{-1}(z_{\ell j}-\mu_j)
\right].$$ Using the Gaussian update, $$\Delta_{\ell j}
=
\rho_j
\left[
\alpha_{\ell j}^2d_j(y_\ell)
+
\frac{J-1}{C_{z,\ell j}}
\right].$$ The update for $q_{\ell j}(\omega_{\ell j})$ is
$$q_{\ell j}(\omega_{\ell j})
=
\operatorname{Gamma}
\left(
\alpha_{\omega,\ell j},
\beta_{\omega,\ell j}
\right),$$ where $$\alpha_{\omega,\ell j}
=
\frac{\nu_0+J+2}{2},$$ and $$\beta_{\omega,\ell j}
=
\frac{1+\Delta_{\ell j}}{2}.$$

For each effect $\ell$, update $$\pi_{\ell j}
=
\frac{\exp\{\mathcal F_{\ell j}\}}
{\sum_{k=1}^J\exp\{\mathcal F_{\ell k}\}},$$ where $$\mathcal F_{\ell j}
=
E_{\ell j}
\left[
\log p_{\nu_0}(s_{\ell j}\mid j)
+
\log p_{\nu_0}(\omega_{\ell j})
+
\log p(z_{\ell j}\mid \omega_{\ell j},j)
+
\log p(\beta_\ell)
+
\log p(y_\ell\mid \beta_\ell,s_{\ell j},z_{\ell j},j)
-
\log q_{\ell j}(s_{\ell j})
-
\log q_{\ell j}(\omega_{\ell j})
-
\log q_{\ell j}(z_{\ell j})
-
\log q_{\ell j}(\beta_\ell)
\right].$$

After updating $\pi_{\ell j}$ and the local factors, update the
posterior mean contribution $$m_\ell
=
\sum_{j=1}^J
\pi_{\ell j}
m_{\beta,\ell j}
\bar s_{\ell j}
\bar u_{\ell j},$$ where $$\bar u_{\ell j}
=
E_{\ell j}[u_j(z_{\ell j})].$$ Let $$g_j=\frac{\bar R_j}{\rho_j}.$$ Then
$g_{jj}=1$, and $$\bar u_{\ell j}
=
g_j
+
\alpha_{\ell j}
\left(
y_\ell-y_{\ell j}g_j
\right).$$

## ELBO for monitoring convergence

Let $$m=\sum_{\ell=1}^L m_\ell.$$ The expected log-likelihood is
$$E_q\log p(x\mid \cdots)
=
-\frac12
\left[
J\log(2\pi)
+
\log|\bar R|
-
J\log N
+
N
\left\{
(x-m)^TK(x-m)
+
\sum_{\ell=1}^L
\left\{
E_q[\theta_\ell^TK\theta_\ell]
-
m_\ell^TKm_\ell
\right\}
\right\}
\right],$$ where $$\theta_\ell
=
\beta_\ell R^{(\ell)}_{\gamma_\ell}.$$ For each effect,
$$E_q[\theta_\ell^TK\theta_\ell]
=
\sum_{j=1}^J
\pi_{\ell j}
E_{\ell j}[\beta_\ell^2]
E_{\ell j}[s_{\ell j}^2]
A_{\ell j}.$$

The full ELBO is $$\mathcal L(q,\nu_0)
=
E_q\log p(x\mid \cdots)
+
\sum_{\ell=1}^L
\sum_{j=1}^J
\pi_{\ell j}
\left[
\log\frac1J
-
\log\pi_{\ell j}
+
\mathcal G_{\ell j}
\right],$$ where $$\mathcal G_{\ell j}
=
E_{\ell j}
\left[
\log p_{\nu_0}(s_{\ell j}\mid j)
+
\log p_{\nu_0}(\omega_{\ell j})
+
\log p(z_{\ell j}\mid \omega_{\ell j},j)
+
\log p(\beta_\ell)
-
\log q_{\ell j}(s_{\ell j})
-
\log q_{\ell j}(\omega_{\ell j})
-
\log q_{\ell j}(z_{\ell j})
-
\log q_{\ell j}(\beta_\ell)
\right].$$ This ELBO can be evaluated after each outer iteration to
monitor convergence.

## Variational empirical Bayes update for $\nu_0$

The empirical Bayes update for $\nu_0$ maximizes the part of the ELBO
depending on $\nu_0$: $$\mathcal M(\nu_0)
=
\sum_{\ell=1}^L
\sum_{j=1}^J
\pi_{\ell j}
\left[
E_{\ell j}\log p_{\nu_0}(s_{\ell j}\mid j)
+
E_{\ell j}\log p_{\nu_0}(\omega_{\ell j})
\right].$$ Since $$s_{\ell j}\mid \nu_0,j
\sim
\operatorname{InvGamma}(a,b_j),$$ with $$a=\frac{\nu_0+2}{2},
\qquad
b_j=\frac{\nu_0\rho_j}{2},$$ we have
$$E_{\ell j}\log p_{\nu_0}(s_{\ell j}\mid j)
=
a\log b_j
-
\log\Gamma(a)
-
(a+1)E_{\ell j}[\log s]
-
b_jE_{\ell j}[1/s].$$ Also, $$\omega_{\ell j}\mid \nu_0
\sim
\operatorname{Gamma}
\left(c,\frac12\right),
\qquad
c=\frac{\nu_0+3}{2},$$ so $$E_{\ell j}\log p_{\nu_0}(\omega_{\ell j})
=
c\log\frac12
-
\log\Gamma(c)
+
(c-1)E_{\ell j}[\log\omega]
-
\frac12E_{\ell j}[\omega].$$ Thus $$\nu_0^{\mathrm{new}}
=
\arg\max_{\nu_0}
\mathcal M(\nu_0).$$ A convenient constrained parameterization is
$$\nu_0=2+\exp(\eta),$$ which enforces $\nu_0>2$.

## Posterior inclusion probabilities

As in SuSiE, the posterior inclusion probability for coordinate $j$ is
$$\operatorname{PIP}_j
=
1-
\prod_{\ell=1}^L
(1-\pi_{\ell j}).$$

## Algorithm summary

The SuSiE-style variational empirical Bayes algorithm is:

1.  Initialize $\nu_0$, $\{\sigma_\ell^2\}_{\ell=1}^L$, and
    $\{\pi_{\ell j}\}$.

2.  Initialize each effect mean $m_\ell$, for example $m_\ell=0$.

3.  Repeat until convergence:

    1.  For each effect $\ell=1,\ldots,L$:

        1.  Form the partial residual $$y_\ell=x-\sum_{k\neq \ell}m_k.$$

        2.  For each candidate $j=1,\ldots,J$, update
            $$q_{\ell j}(z_{\ell j}),
            \qquad
            q_{\ell j}(\beta_\ell),
            \qquad
            q_{\ell j}(s_{\ell j}),
            \qquad
            q_{\ell j}(\omega_{\ell j}).$$

        3.  Update $$\pi_{\ell j}
            \propto
            \exp\{\mathcal F_{\ell j}\}.$$

        4.  Update the effect mean $$m_\ell
            =
            \sum_{j=1}^J
            \pi_{\ell j}
            m_{\beta,\ell j}
            \bar s_{\ell j}
            \bar u_{\ell j}.$$

    2.  Update $\nu_0$ by maximizing $$\mathcal M(\nu_0).$$

    3.  Optionally update the effect-size variances:
        $$\sigma_\ell^{2,\mathrm{new}}
        =
        \sum_{j=1}^J
        \pi_{\ell j}
        E_{\ell j}[\beta_\ell^2].$$

    4.  Evaluate the ELBO $\mathcal L(q,\nu_0)$.

4.  Return $\widehat{\nu}_0$, $\{\pi_{\ell j}\}$, and the PIPs.

"'
