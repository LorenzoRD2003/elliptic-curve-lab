# Fundamental Domain Reduction

Source: [src/elliptic_curves/analytic/fundamental_domain.rs](../../src/elliptic_curves/analytic/fundamental_domain.rs)

This is the iterative modular-reduction loop for `τ ∈ ℍ`. The current
implementation alternates between two educational moves:

- translate by an integer to enter the centered strip
- apply `S(τ) = -1/τ` when `|τ| < 1`

and records every step together with the accumulated modular matrix.

```mermaid
flowchart TB
    A["Input: upper-half-plane point τ and max_steps"] --> B["Initialize current = τ, accumulated_matrix = I, steps = []"]
    B --> C{"current already lies in the standard fundamental domain?"}
    C -->|"yes, and no steps yet"| D["Return report with status = AlreadyReduced"]
    C -->|"yes, after some steps"| E["Return report with status = Reduced"]
    C -->|"no"| F{"steps.len() >= max_steps?"}
    F -->|"yes"| G["Return report with status = StepLimitReached"]
    F -->|"no"| H["Compute next_reduction_step(current)"]

    H --> I{"Real part outside the centered strip?"}
    I -->|"yes"| J["Apply translation matrix T^-shift"]
    I -->|"no"| K{"|τ| < 1?"}
    K -->|"yes"| L["Apply S matrix"]
    K -->|"no"| Z0["Return NumericalComparisonFailed"]

    J --> M["Update current = step_matrix(current)"]
    L --> M
    M --> N["Compose accumulated_matrix = step_matrix ∘ accumulated_matrix"]
    N --> O["Push one reduction step into the history"]
    O --> C
```
