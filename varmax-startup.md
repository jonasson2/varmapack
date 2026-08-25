# VARMAX startup

Consider the VARMAX$(p,q,s)$ model

$$
\tag{1}
x_t =
\sum_{i=1}^p A_i x_{t-i}
+ \varepsilon_t
+ \sum_{i=1}^q B_i\varepsilon_{t-i}
+ \sum_{i=1}^s D_i z_{t-i+1},
\qquad \varepsilon_t\sim\mathcal N(0,\Sigma),
$$

where the exogenous sequence $z_t$ is known. Since this sequence is fixed, it
does not affect the covariance structure, and the ordinary VARMA covariance
results remain applicable. Let

$$
t_0=\max(p,s-1),
\qquad h\geq\max(p,q,s-1).
$$

Let $v=\varepsilon_{t_0-q:h-1}$. Writing $B_0=I$, define the
$1\times(h-t_0+q)$ block matrix of $r\times r$ blocks

$$
H_t=\left[H_{t\ell}\right]_{\ell=t_0-q}^{h-1},
$$

where

$$
H_{t\ell}=
\begin{cases}
B_{t-\ell}, & 0\leq t-\ell\leq q,\\
0,          & \text{otherwise},
\end{cases}
\qquad \ell=t_0-q,\ldots,h-1.
$$

For $t=t_0,\ldots,h-1$, model (1) imposes the constraints

$$
\tag{2}
H_tv =
x_t-\sum_{i=1}^p A_i x_{t-i}
    -\sum_{i=1}^s D_i z_{t-i+1}.
$$

Let $H$ be the matrix with block rows $H_{t_0},\ldots,H_{h-1}$, and let $b$
be the vector obtained by stacking the right-hand sides of equation (2). The
constraints can then be written as $Hv=b$.

Let

$$
V=\operatorname{Var}(v)=I_{h-t_0+q}\otimes\Sigma.
$$

Since $b=Hv$, the vectors $v$ and $b$ are jointly Gaussian with zero means and
$\operatorname{Var}(v)=V$, $\operatorname{Var}(b)=W=HVH^T$, and
$\operatorname{Cov}(v,b)=VH^T$. It follows that

$$
\tag{3}
v\mid b\sim\mathcal N(e_b,R_b),
\qquad e_b=VH^TW^{-1}b,
\qquad R_b=V-VH^TW^{-1}HV.
$$

A realization of $v=\varepsilon_{t_0-q:h-1}$ can therefore be drawn from
equation (3), and it will satisfy the constraints in equation (2). The
simulation can then be carried on with equation (1).
