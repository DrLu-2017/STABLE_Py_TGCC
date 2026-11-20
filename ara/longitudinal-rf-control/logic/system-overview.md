# System overview

## Objective

Describe the conceptual placement of longitudinal RF control features inside `mbtrack2` so that each new feature is implemented as part of a coherent system rather than as an isolated patch.

## System boundary

The relevant system boundary for this artifact includes:

- longitudinal beam dynamics;
- bunch centroid or collective longitudinal observables;
- RF reference quantities and cavity-related state updates;
- feedback/controller objects acting on longitudinal variables;
- perturbation sources such as RF phase noise.

## Conceptual signal chain

The intended signal chain is:

1. beam state evolves turn by turn;
2. a longitudinal observable is measured or estimated;
3. optional controller logic processes the signal;
4. controller output modifies an RF-related actuator or reference;
5. the updated RF action affects the next beam evolution.

## Main extension families

### 1. Phase loop mode 0

This feature targets the global synchronous phase motion associated with mode 0. It should define a measurable phase-like error signal, a controller law, and an explicit actuation point.

### 2. DLLRF

This feature introduces a more detailed low-level RF control representation, potentially including controller dynamics, delay, filtering, and interactions between cavity state and reference signals.

### 3. RF phase noise

This feature introduces a stochastic or prescribed phase perturbation source into the RF chain and studies how it propagates to beam observables.

## Shared architectural principles

- Every feature must expose its measurement path.
- Every feature must expose its actuation or injection path.
- Every feature must define configuration parameters explicitly.
- Every feature must provide at least one isolated validation case.
- Shared infrastructure should be reused only when it preserves physical clarity.

## Non-goals

This top-level artifact does not assume that all longitudinal control and noise features must be merged into one monolithic controller. Separate modules are preferred when they have different physical roles.
