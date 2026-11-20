# Project instructions for GitHub Copilot

## Repository context

This repository contains ongoing work on longitudinal RF control in `mbtrack2`.

Important project fact:
- `rf.py` already contains initial implementations or partial class structures for **DLLRF** and **mode 0 damper**.
- These implementations are not yet fully validated.
- Current work should therefore focus on **understanding**, **reusing**, **validating**, and **minimally extending** the existing code.
- Do **not** start by redesigning these modules from scratch.

## Working rules

1. Always inspect existing classes in `rf.py` before proposing a new controller class.
2. Prefer extension, cleanup, validation hooks, tests, and minimal refactors over replacement.
3. If proposing a new class or a broad refactor, explain explicitly why the existing implementation in `rf.py` is insufficient.
4. Treat `rf phase noise` as a separate future module unless the current task explicitly asks for it.
5. Keep **DLLRF**, **mode 0 damper**, and **RF phase noise** conceptually separate unless a shared interface is clearly justified.
6. Avoid hidden changes to simulation behavior. Make measurement path, controller logic, and actuation path explicit.
7. When modifying architecture, update the ARA trace files as part of the task.

## ARA files to read first

Before making architectural suggestions, read these files if they exist:

- `ara/longitudinal-rf-control/README.md`
- `ara/longitudinal-rf-control/logic/system-overview.md`
- `ara/longitudinal-rf-control/logic/integration-roadmap.md`
- `ara/longitudinal-rf-control/src/index.md`
- relevant submodule files under:
  - `ara/longitudinal-rf-control/phase-loop-mode0/`
  - `ara/longitudinal-rf-control/dllrf/`
  - `ara/longitudinal-rf-control/rf-phase-noise/`

## Default workflow for coding tasks

For tasks related to DLLRF or mode 0 damper:

1. Read `rf.py` first.
2. Identify existing classes, methods, and data flow.
3. Summarize what already exists before suggesting code changes.
4. Identify validation gaps separately from implementation gaps.
5. Propose the smallest safe patch.
6. Only after that, suggest refactors if they materially improve clarity or reuse.

## Response style for repository tasks

When asked to help with DLLRF or mode 0 damper, structure the response as:

1. Existing implementation found in `rf.py`
2. What is missing or unclear
3. Minimal recommended modification path
4. Validation plan
5. Optional future refactor ideas

## Constraints

- Do not duplicate an existing class under a new name without justification.
- Do not create a parallel feedback architecture if the existing one can be reused.
- Do not mix unvalidated DLLRF changes with unrelated RF phase noise work in the same patch.
- Do not replace working code with speculative abstractions.

## Preferred development mindset

This repository is in an incremental research-engineering stage.
The priority is not to invent a fresh software design, but to:

- recover the intent of the current implementation,
- validate it,
- document decisions,
- and evolve it carefully.
