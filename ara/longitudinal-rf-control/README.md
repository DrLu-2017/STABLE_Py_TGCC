# Longitudinal RF Control ARA

This artifact organizes planned and implemented extensions to `mbtrack2` related to longitudinal RF control and perturbation modeling. The current roadmap covers three staged additions:

1. phase loop feedback for mode 0;
2. DLLRF modeling and integration;
3. RF phase noise injection and response analysis.

## Purpose

The purpose of this artifact is to make the design intent, code touchpoints, implementation history, and validation evidence legible to both human collaborators and AI coding agents.

## How to use this artifact

- Read `logic/system-overview.md` first for the conceptual architecture.
- Read `logic/integration-roadmap.md` for staged scope.
- Read `src/index.md` for code navigation.
- Read `trace/master-decisions.yaml` before modifying architecture.
- Enter the relevant submodule directory before working on a specific feature.

## Modules

- `phase-loop-mode0/`: control of bunch-synchronous longitudinal phase motion.
- `dllrf/`: explicit low-level RF loop modeling and associated controller structure.
- `rf-phase-noise/`: noise-source definition, injection path, and beam response validation.

## Policy

Every design change must update:

- one logic file,
- one trace file,
- and one evidence entry if behavior changed.
