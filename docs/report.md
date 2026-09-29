# Report — Práctica 4 (Markdown version)

> Structured English summary of [`original/Practica4.pdf`](original/Practica4.pdf)
> (6 pages, dated Friday 22 February 2019). The PDF remains the authoritative original.
> Figures below were extracted unmodified from the PDF.

## 1. Introduction (page 1)

Evolution Strategies (ES) are a kind of Evolutionary Algorithm mainly characterized by:

- selection of individuals for recombination is unbiased, and selection is a deterministic process;
- they differ from other evolutionary algorithms mainly in the form of the **mutation operator**;
- they are applied mainly to **continuous optimization** problems where individuals are real-valued vectors;
- they were originally created at the Technical University of Berlin in 1964.

The general notation is (μ/ρ +, λ)-ES, where:

- **μ** — population size
- **ρ** — number of parents selected for recombination
- **λ** — number of offspring

*Note:* the PDF text says "has the following notation:" but no formula is rendered after it; the
notation above is the standard one implied by the listed symbols (Inferred).
Source cited in the report: Spanish Wikipedia, *Estrategia evolutiva*.

## 2. Results — "first function" (pages 2–3)

The plots (axis range −20…20, values up to ~800, paraboloid shape) correspond to the **sphere
function** Σ (xᵢ − 2)², i.e. the *second* function of the assignment (see [Discrepancies](#5-discrepancies)).

### (μ + λ)-ES — μ = 100, λ = 30

| Gen | Figure | Author's observation |
|---|---|---|
| 1 | ![](../assets/images/report/sphere-plus-gen01.jpeg) | Initial population is spread within the bounds. |
| 10 | ![](../assets/images/report/sphere-plus-gen10.jpeg) | Population converged to the optimum by generation 4. |
| 15 | ![](../assets/images/report/sphere-plus-gen15.jpeg) | Mutations are now confined to a very small range, unlike the other method, where mutations push children out of the converged point. |

### (μ, λ)-ES — μ = 100, λ = 90

| Gen | Figure | Author's observation |
|---|---|---|
| 1 | ![](../assets/images/report/sphere-comma-gen01.jpeg) | Initial population is spread within the bounds. |
| 10 | ![](../assets/images/report/sphere-comma-gen10.jpeg) | Converged very fast, almost from generation 3. |
| 15 | ![](../assets/images/report/sphere-comma-gen15.jpeg) | Population keeps varying; more individuals move away from the solution due to mutation. |

## 3. Results — "second function" (pages 4–5)

Plots (axis range −4…4, values ±0.4, one valley and one peak) correspond to
**f(x, y) = x·e^(−x² − y²)**.

### (μ + λ)-ES — μ = 100, λ = 30

| Gen | Figure | Author's observation |
|---|---|---|
| 1 | ![](../assets/images/report/xexp-plus-gen01.jpeg) | Initial population is spread within the bounds. |
| 10 | ![](../assets/images/report/xexp-plus-gen10.jpeg) | *(no comment in the report)* |
| 15 | ![](../assets/images/report/xexp-plus-gen15.jpeg) | Mutated children struggle to get in and converge. |

### (μ, λ)-ES — μ = 100, λ = 90

| Gen | Figure | Author's observation |
|---|---|---|
| 1 | ![](../assets/images/report/xexp-comma-gen01.jpeg) | Initial population is spread within the bounds. |
| 10 | ![](../assets/images/report/xexp-comma-gen10.jpeg) | Even after 10 generations some individuals keep mutating outside the valley ("hoyo"). |
| 15 | ![](../assets/images/report/xexp-comma-gen15.jpeg) | No significant changes between generations 5 and 15 (top view). |

## 4. Conclusions (page 6)

To decide which strategy is better, the problem must be defined first. From the experiments:

- **(μ + λ)-ES** evolves more slowly, but once close to the solution individuals barely vary.
- **(μ, λ)-ES** converges extremely fast, but near the solution mutated children move away every
  generation, so some offspring always stay away from the solution.
- **Recommendation:** for (μ, λ)-ES keep λ (children born per generation) at **at most 20% of μ**,
  so fewer solutions diverge from the result near convergence — because mutation (normal-distribution
  based) is always applied; there is no mutation probability.

## 5. Discrepancies

Documented as found; not resolved.

| # | Discrepancy | Sources |
|---|---|---|
| 1 | Function order: the report's "first function" is the sphere function, which is the *second* in the assignment (and vice versa). | Assignment screenshot vs. plot shapes/axes on pages 2–5 |
| 2 | Page 5, generation 10 of (μ, λ)-ES lists λ = 30, while the rest of the section lists λ = 90. Likely a typo (Inferred). | Page 5 |
| 3 | Assignment domain for x·e^(−x²−y²) is [−2, 2]; the code initializes in [−5, 5]. | Assignment vs. `src/*.m` |
| 4 | The sphere-function surface is plotted on the grid `meshgrid(-20:.5:20,-10:.5:20)` while individuals are initialized in [−5, 5]. | Commented block in `src/*.m` |
