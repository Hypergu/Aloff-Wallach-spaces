# Eta-invariant computations for Aloff–Wallach spaces

Symbolic computations supporting the paper

> A. Aleshin, *The Crowley–Nordström invariants of $G_2$-structures on Aloff–Wallach spaces*, https://arxiv.org/abs/2604.04605

The notebook `Computing_eta_inv_AF_spaces.ipynb` carries out the computations of Section 5.2 (*Spectra of deformed Dirac operators*) of the paper, for the Aloff–Wallach space $N_{k,l} = SU(3)/S^1_{k,l}$. The output of the notebook computations supports the Lemma 5.2 from the paper: the spectral flow terms are $J_D = 0$ and $J_B = -2$.

The odd signature operators lie in a one-parameter family

$$ {}^\gamma B^{\lambda,\,3\lambda-1} \;=\; {}^\gamma\widetilde{B} + \mu B_0, \qquad \mu = 3\lambda - 1, \quad \lambda \in [\tfrac13, \tfrac12], $$

connecting the reductive operator $\widetilde{B} = B^{1/3,0}$ at $\mu = 0$ to the Levi-Civita odd signature operator $B^{1/2,1/2}$ at $\mu = \tfrac{1}{2}$. So $B_0$ is the *difference* between the two operators. There are two stages of computations:

1. The spectrum of $B_0$ on $S \otimes S$ bounds how far the deformation can move an eigenvalue. Comparing $\lambda_{\max}\big(\tfrac{1}{2} B_0\big)$ against the lower bound for $({}^\gamma\widetilde{B})^2$ shows that no eigenvalue can cross zero along the path except in the **trivial** representation $\gamma_0$.

2. For $\gamma_0$, both endpoint operators are restricted to the 10-dimensional subspace $\mathrm{Hom}_H\big(V^{(0,0)}, S\otimes S\big)$ and their spectra computed directly, giving $\eta$ and $h$ for each and hence $J_B$.

Everything is done in symbolic arithmetic in $k$ and $l$ — there is no numerics and no floating point computations.

## Requirements

- Python 3
- [SymPy](https://www.sympy.org)

## What the notebook does

1. **Basis of $\mathfrak{su}(3)$ adapted to $N_{k,l}$.** We set up $e_1,\dots,e_7$, with $e_4$ spanning $\mathfrak{s}^\perp$, and print the commutator table. This reproduces the commutator relations recorded in Section 4.1 of the paper.

2. **Spinor model.** We identify the spinor space $S$ with $\Lambda^*(\mathbb{C}^3)$ in the basis $(1,\ f_1,\ f_2,\ f_3,\ f_1\wedge f_2,\ f_2\wedge f_3,\ f_3\wedge f_1,\ f_1\wedge f_2\wedge f_3)$, written $(s_1,\dots,s_8)$. The operators $\varepsilon_i$ (left multiplication by $f_i$) and $\iota_i$ (contraction by the dual vector) are given as explicit $8\times 8$ matrices, from which the Clifford multiplications $c_1,\dots,c_7$ are built.

3. **The deformation term $B_0$.** We represent $B_0 = \sum_i c_i\big(\frac{1}{3}\widetilde{\mathrm{ad}}_i + \hat{\widetilde{\mathrm{ad}}}_i\big)$ as a $64\times 64$ matrix on $S\otimes S$, as `B_left + B_right`.

4. **Spectrum of $B_0$.** We split $B_0$ into blocks by connected components of its nonzero pattern and compute each block's eigenvalues. The largest eigenvalue has absolute value $8\sqrt{2}$, i.e. $2\sqrt{2}$ after applying the paper's normalisation (see *Normalisation* below), so $\lambda_{\max}\big(\tfrac12 B_0\big)^{2} = 2$ — the constant the paper compares against the lower bound for $({}^\gamma\widetilde{B})^{2}$ to rule out every representation but $\gamma_0$.

5. **The trivial representation: $\widetilde{B} = {}^{\gamma_0}B^{1/3,0}$.** We restrict `B_left` to the 10-dimensional subspace $\mathrm{Hom}_H\!\big(V^{(0,0)}, S\otimes S\big)$, spanned by

   $$s_1\otimes s_1,\; s_1\otimes s_8,\; s_2\otimes s_6,\; s_3\otimes s_7,\; s_4\otimes s_5,\; s_5\otimes s_4,\; s_6\otimes s_2,\; s_7\otimes s_3,\; s_8\otimes s_1,\; s_8\otimes s_8,$$

   which is what the index list `H_inv` encodes. Every eigenvalue is nonzero, giving

   $$\eta\big({}^{\gamma_0}B^{1/3,0}\big) = 2\,\mathrm{sign}(k) + 2\,\mathrm{sign}(l) - 2\,\mathrm{sign}(k+l) = 2, \qquad h\big({}^{\gamma_0}B^{1/3,0}\big) = 0$$

   under the assumptions below.

6. **The trivial representation: ${}^{\gamma_0}B^{1/2,1/2}$.** We form `S = B_left + B0/2` and restrict it to the same $H$-invariant subspace. Some blocks have eigenvalues SymPy cannot express in closed form, but only their *signs* are needed for our applications: each block is self-adjoint, so its eigenvalues are real and the sign counts follow from Descartes' rule applied to the characteristic polynomials. The $10\times10$ block is handled by factoring its characteristic polynomial as

   $$x^{2}\,(x^{2}-96)\,(x^{2}-48)\Big(x^{4}-144x^{2}+96\sqrt{6}\,\tfrac{kl(k+l)}{(k^{2}+kl+l^{2})u}\,x+3072\Big), \qquad u=\sqrt{k^{2}+kl+l^{2}},$$

   whose quartic factor has two positive and two negative roots, giving $\eta\big({}^{\gamma_0}B^{1/2,1/2}\big) = 0$ and $h\big({}^{\gamma_0}B^{1/2,1/2}\big) = 2$.

7. **Conclusion.** Combining the two,

   $$\eta\big({}^{\gamma_0}B^{1/2,1/2}\big) - (\eta+h)\big({}^{\gamma_0}B^{1/3,0}\big) = 0 - 2 = -2,$$

   which is the $J_B = -2$ of Lemma 5.2.

## Assumptions

- $k, l$ are declared **positive integers**. The paper's results cover general pairs $(k,l)$; the remaining cases follow from the sign and symmetry arguments for the Aloff–Wallach spaces.
- $k \neq l$ is used in step 6: one block's sign is determined by $\mathrm{sign}(-k^{2}+2kl-l^{2}) = \mathrm{sign}\big(-(k-l)^{2}\big)$, which is negative precisely when $k \neq l$. This condition is automatically satisfied in the generic case.
- The notebook treats the *generic* case. The exceptional spaces $N_{1,1}$ and $N_{1,0}$ have a larger $H$-invariant subspace of $S\otimes S$ and are handled separately in the paper's appendix. However, they can be computed from the notebook's code — one just needs to use the basis for the $H$-invariant subspace in each exceptional case.

## Normalisation

- The $\widetilde{\mathrm{ad}}_i$ are **not** divided by 4, because doing so makes SymPy use floating point. Every eigenvalue printed for $B_0$ is therefore 4× the paper's value: the maximum $8\sqrt{2}$ corresponds to $2\sqrt{2}$ there.

## License

MIT — see `LICENSE`.
