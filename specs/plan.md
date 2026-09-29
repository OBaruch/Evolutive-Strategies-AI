# Plan

> Two parts: **(A)** the original development plan, reconstructed from the artifacts (the order is
> Inferred), and **(B)** the plan executed later to reorganize the repository.

## A. Original implementation plan (reconstructed, 2019)

| Step | Activity | Artifact |
|---|---|---|
| A1 | Read the assignment and study ES theory (μ, ρ, λ notation). | Report §Introduction |
| A2 | Implement the (μ + λ)-ES: init → recombine → mutate → evaluate → sort → plot. | `src/UmasLambda.m` |
| A3 | Copy it and change evaluation/selection to obtain the (μ, λ)-ES. | `src/UcomaLambda.m` |
| A4 | Run both scripts on the sphere function, capture generations 1, 10 and 15. | Report pp. 2–3 |
| A5 | Switch the active function to x·e^(−x²−y²) and repeat. | Report pp. 4–5 (the committed code keeps this function active) |
| A6 | Write conclusions and export the report to PDF. | Report p. 6, `docs/original/Practica4.pdf` |
| A7 | (2021) Upload scripts, report and MIT license to GitHub. | Git history |

## B. Repository reorganization plan (executed)

Principle: **modernize the repository, not the project.**

| Step | Task | Done |
|---|---|---|
| B1 | Inventory all files; read code, PDF text and inspect every PDF page visually. | ✅ |
| B2 | Classify the project with evidence (Academic / Coursework). | ✅ |
| B3 | Move code to `src/` and the report to `docs/original/` with `git mv` (content unchanged). | ✅ |
| B4 | Extract report figures and assignment screenshots into `assets/images/report/`. | ✅ |
| B5 | Write `README.md` and `docs/` (context, assignment, report, code overview, improvements). | ✅ |
| B6 | Write `specs/` (intent, spec, plan) and `AGENTS.md`. | ✅ |
| B7 | Add a minimal MATLAB `.gitignore`. | ✅ |
| B8 | Verify `src/` blobs are identical to the original commit. | ✅ |

### Verification

```bash
# Must print nothing: the original blobs and the moved blobs are identical
for f in UmasLambda.m UcomaLambda.m; do
  [ "$(git rev-parse 4d6f354:$f)" = "$(git rev-parse HEAD:src/$f)" ] || echo "CHANGED: $f"
done
```

## C. Out of scope

- Any change to the MATLAB code (fixes listed in `docs/possible-improvements.md` are not applied).
- CI/CD, tests, containers, package managers or other tooling.
