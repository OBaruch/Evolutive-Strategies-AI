# AGENTS.md

Guidance for anyone — human contributor or AI coding agent — working on this repository.

## Nature of the repository

A **historical, preserved** university project (2019). Its value is the original implementation.
Read [`README.md`](README.md) and [`specs/`](specs) before changing anything.

## Hard rules

1. **Never modify files in `src/`.** No fixes, formatting, renames, line-ending changes or comment edits.
2. **Never modify files in `docs/original/`.**
3. Do not add build systems, CI/CD, containers, linters or test frameworks.
4. Improvements go to [`docs/possible-improvements.md`](docs/possible-improvements.md) as text only.
   A future rewrite, if any, must live in a separate directory and be clearly labeled as new work.
5. In documentation, label claims as **Confirmed**, **Inferred** or **Unknown**; do not invent context.

## Workflow

`specs/intent.md` → `specs/spec.md` → `specs/plan.md` → change → verify that `src/` is unchanged:

```bash
# Must print nothing: the original blobs and the moved blobs are identical
for f in UmasLambda.m UcomaLambda.m; do
  [ "$(git rev-parse 4d6f354:$f)" = "$(git rev-parse HEAD:src/$f)" ] || echo "CHANGED: $f"
done
```
