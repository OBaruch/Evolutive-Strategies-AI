# Specification

> Reverse-engineered from [`../src`](../src) and [`../docs/original/Practica4.pdf`](../docs/original/Practica4.pdf).
> Describes the system **as built**. Where the implementation deviates from the assignment or from
> textbook ES, the deviation is recorded, not corrected.

## 1. Scope

Two MATLAB scripts that minimize a 2-D function with an Evolution Strategy and animate the
population.

| ID | Component | File |
|---|---|---|
| C1 | (μ + λ)-ES | `src/UmasLambda.m` |
| C2 | (μ, λ)-ES | `src/UcomaLambda.m` |

## 2. Functional requirements

| ID | Requirement | Source | Status in code |
|---|---|---|---|
| FR-1 | Minimize f(x, y) = x·e^(−x²−y²). | Assignment | Implemented (active block) |
| FR-2 | Minimize f(x) = Σ(xᵢ − 2)², d = 2. | Assignment | Implemented (commented block, toggled manually) |
| FR-3 | Provide a (μ + λ)-ES. | Assignment | C1 |
| FR-4 | Provide a (μ, λ)-ES. | Assignment | C2 (selection deviates, see §5) |
| FR-5 | Initialize μ individuals uniformly within bounds with random σ ∈ [0,1). | Code | C1, C2 |
| FR-6 | Create offspring by intermediate recombination of 2 distinct random parents (genes and σ). | Code | C1, C2 |
| FR-7 | Mutate offspring with Gaussian noise N(0, σ). | Code | C1, C2 |
| FR-8 | Rank individuals ascending by fitness (minimization). | Code | C1, C2 |
| FR-9 | Animate population over the function's `surfc` plot each generation. | Code / report | C1, C2 |
| FR-10 | Deliver a PDF report with figures and conclusions. | Assignment | `docs/original/Practica4.pdf` |

## 3. Parameters

| Name | Meaning | C1 | C2 |
|---|---|---|---|
| `G` | generations | 15 | 15 |
| `mu` | population size | 100 | 100 |
| `lambda` | offspring per generation | 30 | 90 |
| `xl`, `xu` | lower / upper bounds | [−5, −5], [5, 5] | same |
| `D` | genes (dimensions) | 2 | 2 |

## 4. Interfaces

- **Input:** none; edit the script. Function switched by comment toggling.
- **Output:** a MATLAB figure animation (0.5 s per generation) and σ printed to the Command Window.
- **Dependencies:** MATLAB + Statistics and Machine Learning Toolbox (`normrnd`).

## 5. Known deviations (accepted as historical)

| ID | Deviation |
|---|---|
| D-1 | Domain [−5, 5] instead of [−2, 2] for FR-1; no bound enforcement. |
| D-2 | Offspring loop starts at index μ (overwrites last parent, λ + 1 children). |
| D-3 | All children are mutated with σ of column μ + λ; σ is not self-adapted. |
| D-4 | C2 selection mixes the λ top-ranked individuals with previous columns instead of selecting only among offspring. |
| D-5 | Two code files, while the assignment asked for one. |

Details: [`../docs/code-overview.md`](../docs/code-overview.md).

## 6. Acceptance (as reported, 2019)

| Scenario | Reported outcome |
|---|---|
| Sphere, (μ + λ), μ=100, λ=30 | Converged by generation 4; stable afterwards. |
| Sphere, (μ, λ), μ=100, λ=90 | Converged ~generation 3; offspring keep drifting. |
| x·e^(−x²−y²), (μ + λ) | Population concentrates in the valley; children struggle to enter. |
| x·e^(−x²−y²), (μ, λ) | Some individuals keep mutating out of the valley; no significant change after gen 5. |

## 7. Repository constraints (current)

- Files under `src/` MUST remain byte-for-byte identical to the original upload.
- Original documents under `docs/original/` MUST NOT be modified.
