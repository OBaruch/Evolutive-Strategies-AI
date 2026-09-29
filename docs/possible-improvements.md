# Possible Improvements

> **None of these improvements have been applied.** The code in [`../src`](../src) is intentionally
> kept exactly as originally written to preserve the historical implementation. This list is only
> a reference for a hypothetical future rewrite, which should live in a separate location.

## Algorithm correctness

| Observation | Where | Possible improvement |
|---|---|---|
| Offspring loop `for t=mu:mu+lambda` overwrites the last parent and creates λ + 1 children. | both | Use `mu+1:mu+lambda`. |
| Mutation uses `sigma(:,mu+lambda)` for every child. | both | Use the child's own `sigma(:,t)`. |
| σ is never mutated (no self-adaptation). | both | Apply log-normal self-adaptation, e.g. `σ' = σ·exp(τ′·N(0,1) + τ·Nᵢ(0,1))`, or the 1/5 success rule. |
| (μ, λ) selection ranks with stale offspring fitness and mixes previous parents into the next generation. | `UcomaLambda.m` | Select the μ best exclusively from the λ offspring (requires λ ≥ μ). |
| No bound handling; the assignment domain for function 1 is [−2, 2] but bounds are [−5, 5]. | both | Match the domain and clamp/reflect offspring. |
| Inner loop reuses the outer loop variable `t`. | both | Use a distinct index (MATLAB still runs G outer iterations, but it harms readability). |

## Code quality

- Extract the ES into a function `es(f, mu, lambda, G, xl, xu, mode)` with `mode ∈ {'plus','comma'}`
  to remove duplication and satisfy the "single code file" deliverable.
- Pass the objective function as a parameter instead of commenting blocks in/out.
- Preallocate `fit`, `xdos`, `sig`; add missing semicolons.
- Replace `normrnd` with `randn .* sigma` to drop the toolbox dependency.
- Seed the RNG (`rng(seed)`) for reproducible experiments.
- Use professional wording in comments.

## Experimentation

- Log best/mean fitness per generation and plot convergence curves instead of relying on screenshots.
- Run multiple seeds and report statistics to compare strategies objectively.
- Report the best solution found vs. the known optimum.
