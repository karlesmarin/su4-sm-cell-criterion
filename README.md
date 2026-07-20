# ⚛️ Three Gates to a Quark Generation

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21432628-1B6F8C?logo=doi&logoColor=white)](https://doi.org/10.5281/zenodo.21432628)
[![License](https://img.shields.io/badge/License-Apache_2.0-B5530F)](LICENSE)
[![Verified](https://img.shields.io/badge/SageMath-exact_to_dim_900-2C2C2C)](https://www.sagemath.org/)
[![Language](https://img.shields.io/badge/paper-EN_%2B_ES-1B6F8C)](.)

**📄 Paper (EN + ES) & all verification scripts on Zenodo → https://doi.org/10.5281/zenodo.21432628**

> ### 📚 Part **II** of a series
> - **Part I — *Anomaly- and Tadpole-Compatible Fermion Completion of 6D SU(4) GHU***
>   → [github.com/karlesmarin/ghu-su4-completion](https://github.com/karlesmarin/ghu-su4-completion) · [Zenodo 10.5281/zenodo.21432625](https://doi.org/10.5281/zenodo.21432625)
> - **Part II — *Three Gates to a Quark Generation*** (this repo): the exact criterion for *which* $SU(4)$ representations contain the SM quark cell.
> - **Part III — *A Centre-Charge Selection Rule for the Wilson-Line Potential***
>   → [github.com/karlesmarin/centre-parity-selection](https://github.com/karlesmarin/centre-parity-selection) · [Zenodo 10.5281/zenodo.21438226](https://doi.org/10.5281/zenodo.21438226)
> - **Part IV — *Schur Functions at $(1,-1,t,t^{-1})$***
>   → [github.com/karlesmarin/schur-nonidentity-o4](https://github.com/karlesmarin/schur-nonidentity-o4) · [Zenodo 10.5281/zenodo.21463000](https://doi.org/10.5281/zenodo.21463000)

**An exact criterion for which $SU(4)$ representations contain the Standard Model.** Part I embedded the
Standard-Model quarks in the dimension-60 representation $(3,\mathbf{60})$ of $SU(4)$ on $T^2/\mathbb{Z}_2$
and found, by an exhaustive scan, that $(3,\mathbf{60})$ is the *smallest* representation able to do so.
This Part II turns that scanned fact into a law.

## 🎯 The main result

An irreducible representation with Dynkin labels $(a,b,c)$ contains the Standard-Model quark cell
$\{Q(\mathbf{2},\tfrac16),\,u(\mathbf{1},\tfrac23),\,d(\mathbf{1},-\tfrac13)\}$ in its $T^2/\mathbb{Z}_2$
chiral zero modes **if and only if**

```
(a + 2b + 3c) odd    ∧    b ≥ 1    ∧    a + b + c ≥ 3
```

— an **arithmetic** condition (the $\mathbb{Z}_4$ centre charge), a **geometric** one (the middle node of
$SU(4)$ excited), and a **size** one (enough room to close the cell). Up to dimension 400 only **eight**
representations qualify, of dimensions $60,84,140,140,216,224,280,360$; the smallest, $(0,2,1)=\mathbf{60}$,
is exactly the $(3,\mathbf{60})$ of Part I — recovered and *explained*.

## 🧩 How it is derived

Delete the middle node of the $SU(4)$ Dynkin diagram $\circ\!-\!\circ\!-\!\circ$ → a Levi subgroup
$SU(2)_L\times SU(2)_R\times U(1)_m$. The right-handed colour-singlet zero modes form $(b{+}1)$ copies of the
Clebsch–Gordan tower $[a/2]\otimes[c/2]$; the orbifold chiral projection keeps half of each, giving the
**closed-form zero-mode count**

```
N(a,b,c) = (b+1)(a+c+1)/2
```

**Honest credit.** The centre gate is classical (Hucks, Saller, Tong, Slansky); the branching is
Littlewood–Richardson; the midpoint property ($\tfrac16$ = mean of $\tfrac23,-\tfrac13$) reflects the
Pati–Salam hypercharge relation. What is ours is the *orbifold chiral projection* that yields the closed
count, the packaged three-gate criterion, and the reading of the $\pm\tfrac12$ cell as the invariant.

## 📈 Why $SU(4)$? The critical rank

The cell imposes three charge constraints on a hypercharge functional with $\mathrm{rank}(G)-1$ free
parameters, so it is a **rigid** obstruction exactly at $\mathrm{rank}=3$ — $SU(4)$ is the unique marginal
case. For $SU(5),SU(6)$ (rank $\ge 4$) a generic representation admits (73% / 85% to dim 200, across all
centre classes). The same counting explains why a *pinned* hypercharge in higher-rank models (Komori–Maru
$SU(7)$) leaves quarks with integer charges and forces them onto a brane.

## ✅ Reproducibility — every count is exact (SageMath), each theorem has its script

There is no Lean certificate here: the statements are ∀-over-all-representations, not a finite `decide`, so
they are machine-checked in **exact rational arithmetic** over $SU(n)$ weight systems. Each script below is
mapped to the theorem it verifies.

**The cell criterion (Theorem 1, §5).**
- `cell_theorem.sage` — the exact $\det=0$ reformulation of "contains the SM cell" (a rank condition on the
  weight diagram), and its match to the direct construction.
- `cell_unify.sage` — that (given odd centre) *admits* $\iff$ the RH colour-singlet lattice is 2-dimensional.
- `harden.sage` — the full criterion verified against the direct construction on **all reps to dimension 900**.
- `audit.sage` — the authoritative audit: the catalog of admitting reps to dim 400 ($60,84,140,140,216,224,280,360$),
  criterion-vs-construction (0 mismatches), $N=(b{+}1)(a{+}c{+}1)/2$ (0 mismatches), extent $=12b$, and the
  physical hypercharges $\tfrac16,\tfrac23,-\tfrac13$.

**The Levi derivation & the closed count $N$ (§4).**
- `branch_derive.sage` — deletes the middle node, branches to $SU(2)_L\times SU(2)_R\times U(1)_m$, and derives
  the three sub-claims ($b{+}1$ levels, $SU(2)_R$ spin $\tfrac{a+c}{2}$, RH parity keeps half).
- `f1_blindaje.sage` — the corrected derivation: the singlet sector is the Clebsch–Gordan **tower**
  $[a/2]\otimes[c/2]$, whose extent gives $N=(b{+}1)(a{+}c{+}1)/2$.
- `lemma_b2.sage` — the middle-node coordinate $m=2q_8+q_{15}$ is $\alpha_2$-conjugate; the singlet extent is $12b$.
- `lemma_c.sage` — the closed count $N=(b{+}1)(a{+}c{+}1)/2$.
- `doublet_coop.sage` — every admitting rep has a **symmetric** cell $Q=\tfrac12(u+d)$ (so $Y(Q)=\tfrac16$ is automatic).
- `interleave.sage` — the LH doublets interleave the RH singlets at half a middle-node cell (why the cell closes).

**Node-independence (Lemma 2, §6).**
- `node_indep_types.sage` — the Weyl-orbit lemma across Lie types: $A,D$ node-independent; $B,C,G,F$ split by
  long/short root (e.g. $F_4\to[21,21,15,15]$).
- `bonus_lemma.sage` — node-independence for $SU(5),SU(6)$; the product form $(a{+}1)(b{+}1)(c{+}1)$ is special to $SU(4)$.
- `c2_formula.sage` — the general $SU(n)$ node-invariant is a Littlewood–Richardson sum (no clean product for $n\ge5$).

**Why $SU(4)$? The critical rank (§7) and the dimension knob (§8).**
- `sun_cell2.sage` — the $SU(4)/SU(5)/SU(6)$ rank scans: the rigid gate dissolves above rank 3 (73% / 85% admit).
- `c4_center.sage` — the SM cell across the rank-3 simple groups $A_3,B_3,C_3$ (free hypercharge).
- `extrapolate2.sage` — the count factorizes as (Levi branching)$\times$(orbifold projection); the $n_3{+}n_4$
  locking (why 5D and 6D counts coincide).

**Figures.**
- `make_cell_figs.py` — regenerates the three data figures (admissibility landscape, real singlet lattice, rank dissolution).

## 📜 License & citation

Released under [Apache 2.0](LICENSE). Please cite the Zenodo record
**[10.5281/zenodo.21432628](https://doi.org/10.5281/zenodo.21432628)**. Author: **Carles Marín**
(independent researcher). The exact computations were carried out and cross-checked with Claude (Anthropic)
as an AI research assistant against a common machine-verifiable ground truth (exact rational arithmetic on
$SU(n)$ weight systems).
