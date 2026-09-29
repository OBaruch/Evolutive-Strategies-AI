# Project Context

## Classification

**Academic / University Project — Coursework / Assignment** (Confirmed)

## Evidence

| Source | Evidence | Status |
|---|---|---|
| `docs/original/Practica4.pdf` — header of every page | Universidad de Guadalajara / CUCEI logo | Confirmed |
| `docs/original/Practica4.pdf` — page 1 | Author "Omar Baruch Morón López", course "Sistemas Inteligentes II", "Practica 4" | Confirmed |
| `docs/original/Practica4.pdf` — page 1 | Embedded screenshot of the assignment statement ("Práctica 4") | Confirmed |
| `docs/original/Practica4.pdf` — metadata | Author "Omar Baruch", created with Microsoft Word for Office 365 on 2019-02-22 | Confirmed |
| `src/*.m` — header comment | `%% OMAR BARUCH MORON lOPEZ` followed by a student ID number | Confirmed |
| Git history | Files uploaded to GitHub on 2021-02-20 ("Add files via upload") | Confirmed |

## Timeline

| Date | Event |
|---|---|
| 2019-02-22 | Report dated and PDF generated (Word export) |
| 2021-02-20 | Scripts, report and MIT license uploaded to GitHub |
| Later | Repository reorganized and documented (this documentation); source code untouched |

The exact date the scripts were written is **Unknown**; they are assumed (Inferred) to be from
February 2019 because the report shows their output.

## Objective

Implement and compare two Evolution Strategy variants — (μ + λ)-ES and (μ, λ)-ES — for finding
the global minimum of two benchmark functions, and answer which one is better and why.
See [assignment.md](assignment.md).

## Scope

- Two standalone MATLAB scripts, one per strategy.
- Two objective functions, switched by commenting code in/out.
- Visual (animated plot) evaluation; no numeric logging, metrics or tests.
- A 6-page PDF report with screenshots and a written conclusion.

## Learning Context (Inferred)

The course (*Intelligent Systems II*) appears to cover evolutionary computation: the report
introduces Evolution Strategies, the (μ/ρ +, λ) notation, and cites the Spanish Wikipedia article
on *Estrategia evolutiva* as its source. The practice number (4) suggests earlier practices covered
other topics, but the repository does not contain them.

## Unknown

- The instructor, semester/term identifier and grading criteria.
- The exact MATLAB version used.
- Whether a separate code submission existed in a different form (the assignment asks for a single
  code file, but two scripts are present — one per strategy).
