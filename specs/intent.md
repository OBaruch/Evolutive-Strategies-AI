# Intent

> Reconstructed retroactively from the existing repository (assignment, report and code).
> It records *why* the project existed; it is not a plan for new work.
> Evidence levels: **Confirmed** (stated in files) · **Inferred** (deduced) · **Unknown**.

## Why

Complete *Práctica 4* of **Sistemas Inteligentes II** (Universidad de Guadalajara, CUCEI, February
2019): learn how Evolution Strategies work by implementing them, and compare the two classic
replacement schemes experimentally. — Confirmed

## Problem

Find the global minimum of two continuous 2-D functions with an evolutionary method instead of a
gradient-based one: — Confirmed

- f(x, y) = x·e^(−x² − y²), x, y ∈ [−2, 2]
- f(**x**) = Σ (xᵢ − 2)², d = 2

## Desired outcome

1. A working implementation of **(μ + λ)-ES** and **(μ, λ)-ES**. — Confirmed
2. Visual evidence of how each population evolves. — Confirmed
3. A reasoned answer to *"Which evolution strategy is better, and why?"* — Confirmed

## Users / audience

- The course instructor evaluating the practice (Inferred).
- Today: readers of a technical portfolio who want to understand the original work.

## Success criteria (as evidenced by the report)

- The population visibly converges toward the minimum of each function.
- Behavioral differences between both strategies are observed and explained.
- A PDF report with figures and conclusions is delivered.

## Non-goals

- General-purpose optimization library, performance tuning, or statistical benchmarking (Inferred —
  none present).
- Modernizing the original code (explicit constraint of the later repository reorganization).
