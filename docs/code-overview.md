# Code Overview

> Describes the original scripts in [`../src`](../src) **as they are**. Nothing here has been
> applied to the code. Line references are to the unmodified files.

Both scripts are standalone MATLAB scripts (no functions, no shared files) and share ~80% of their
code. They differ only in `lambda` and in the **evaluation/selection** step.

## Shared structure

| Section (original `%%` heading) | Purpose |
|---|---|
| `clear all / close all / clc` | Reset workspace, figures and console. |
| `%% OMAR BARUCH MORON lOPEZ` | Author header (name + student ID). |
| `%% Funcion a analizar` | Objective function `f` as an anonymous function, a `meshgrid` + `z` for the surface, and `axis`. The sphere function block is commented out; `x.*exp(-x.^2-y.^2)` is active. |
| `%% VARIABLES QUE SE PUEDEN MODIFICAR` | Tunable parameters: `G` (generations = 15), `mu` (population = 100), `lambda` (offspring). |
| (unnamed) | Bounds `xl = [-5; -5]`, `xu = [5; 5]`; `D = 2` genes. |
| (unnamed) | Preallocation of `x` (D × (μ+λ)), `sigma` (D × (μ+λ)), `fitness` (1 × (μ+λ)). |
| `%% Inizializar poblacion` | For i = 1…μ: uniform random position in bounds; `sigma = rand(D,1)`. |
| `%% Repetir Generaciones` | Main loop over `t = 1:G`. |

### Main loop, common steps

1. **Offspring creation** — inner loop `for t=mu:mu+lambda` (reuses the name `t`):
   - pick two different random parents `r1 ≠ r2` from 1…μ;
   - intermediate recombination: child genes and σ are the average of both parents;
   - mutation: `r = normrnd(0, sigma(:,mu+lambda))`, `x(:,t) = x(:,t) + r`.
2. **Evaluation and selection** — differs per script (below).
3. **Visualization** — `cla`; for each individual draw a black `*` with `plot3`, plus a red `o`
   for columns after μ; draw `surfc(xfuncion, yfuncion, z)`; `pause(.5)`.

Representation per individual: position `x(:,i)` ∈ ℝ² and a per-gene step size `sigma(:,i)` ∈ ℝ².
Fitness is the raw function value; sorting is ascending, so the goal is **minimization**.

## `UmasLambda.m` — (μ + λ)-ES

- `lambda = 30`.
- Original comment: *"This algorithm reproduces lambda children, mixes them with the parents and
  removes the worst lambda individuals."*
- Evaluation: fitness of **all** μ + λ columns.
- Selection: `sort(fitness)` ascending and reorder `fitness`, `x`, `sigma`. The first μ columns
  (the best of parents ∪ offspring) act as the parents of the next generation; the last columns
  are overwritten by new offspring. This matches the (μ + λ) scheme.

## `UcomaLambda.m` — (μ, λ)-ES

- `lambda = 90`.
- Original comment: *"This algorithm makes lambda children be born and replaces the worst lambda
  parents with them."*
- Evaluation: copies `fitness2 = fitness`, `x2 = x`, `sigma2 = sigma`; evaluates parents
  (1…μ) into `fitness2` and offspring (μ…μ+λ) into `fitness`.
- Selection: sorts `fitness2`, takes the first λ entries into `fit`, `xdos`, `sig`, then builds
  the new population as `[best λ, previous columns 1…μ]`.
- Observed effect (Inferred from reading the code): the next generation's parents are the λ
  top-ranked individuals followed by the first μ − λ columns of the previous array. The ranking
  uses `fitness2`, whose offspring entries hold values from the previous generation (zeros in the
  first one). This differs from the textbook (μ, λ) scheme, where parents are selected only among
  the λ offspring; it is preserved as written.

## Execution flow

```
init f, grid, parameters
init population (μ random individuals)
repeat G times:
    create offspring (recombine 2 parents + Gaussian mutation)
    evaluate
    select  (script-specific)
    redraw population + surface, pause 0.5 s
```

## Dependencies observed

| Call | Origin |
|---|---|
| `normrnd` | Statistics and Machine Learning Toolbox |
| `rand`, `randi`, `sort`, `zeros`, `meshgrid`, `exp` | Core MATLAB |
| `plot3`, `surfc`, `axis`, `hold`, `cla`, `pause` | Core MATLAB graphics |

## Behavioral notes (documented, not changed)

- The inner offspring loop starts at column `mu`, so it overwrites the last parent and produces
  λ + 1 children.
- Mutation uses the step size of column `mu+lambda` for every child, not the child's own σ. In the
  first generation that column is zero until the last child is created.
- σ is recombined but never mutated (no self-adaptation).
- Offspring are not clamped to `[xl, xu]`.
- `sigma = sigma(:,ind)` / `sigma2 = sigma(:,ind)` have no trailing `;`, so σ is printed each generation.
- Files use CRLF line endings.

See [possible-improvements.md](possible-improvements.md) for how these could be addressed in a
future, separate implementation.
