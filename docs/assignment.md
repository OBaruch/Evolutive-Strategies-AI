# Assignment — Práctica 4

> Source: screenshots of the assignment statement embedded in the original report
> ([page 1](../assets/images/report/assignment-statement.png) and
> [page 6](../assets/images/report/assignment-deliverables.png) of
> [`original/Practica4.pdf`](original/Practica4.pdf)). Translated from Spanish.

**Course:** Sistemas Inteligentes II — Universidad de Guadalajara, CUCEI

## Statement

Write a computer program that finds the **global minimum** of the following functions using the
Evolution Strategies **(μ + λ)-ES** and **(μ, λ)-ES**:

1. f(x, y) = x · e^(−x² − y²), with x, y ∈ [−2, 2]
2. f(**x**) = Σᵢ₌₁ᵈ (xᵢ − 2)², with d = 2

## Question to answer

> Which evolution strategy is better? Why?

## Deliverables

- A report in **PDF** format including graphs, tables, etc. showing the results obtained.
- The computer program. Generate **only one file** for the code (`*.m`, `*.c`, `*.cpp`, `*.py`, etc.).
- Do not generate header files.

## Reference: expected optima (Inferred, standard results — not stated in the assignment)

| Function | Global minimum |
|---|---|
| x · e^(−x² − y²) | f ≈ −0.4289 at (x, y) = (−1/√2, 0) ≈ (−0.7071, 0) |
| Σ (xᵢ − 2)² | f = 0 at (2, 2) |

## Compliance notes (Observed)

| Requirement | What the repository contains |
|---|---|
| Both strategies | Yes — one script per strategy (`UmasLambda.m`, `UcomaLambda.m`) |
| Both functions | Yes — switched by commenting/uncommenting code |
| Domain x, y ∈ [−2, 2] for function 1 | Code initializes in [−5, 5] (`xl`, `xu`) and does not clamp |
| Single code file | Two files are present |
| PDF report | Yes — `original/Practica4.pdf` |
