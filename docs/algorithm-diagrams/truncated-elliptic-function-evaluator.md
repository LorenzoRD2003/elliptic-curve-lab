# Truncated Elliptic-Function Evaluator

Source: [src/elliptic_curves/analytic/elliptic_functions/evaluator.rs](../../src/elliptic_curves/analytic/elliptic_functions/evaluator.rs)

This is the shared reduction-and-summation pipeline behind the current
truncated analytic elliptic-function evaluations.

```mermaid
flowchart TB
    A["Input: lattice Λ, point z, truncation, singular_term, lattice_term, builder"] --> B["Reduce z modulo Λ to a canonical fundamental coordinate"]
    B --> C["Recover canonical_z from that coordinate"]
    C --> D{"canonical_z is too close to 0?"}
    D -->|"yes"| Z0["Return PointTooCloseToLatticePoint"]
    D -->|"no"| E["Enumerate all nonzero lattice points ω in the square truncation box"]
    E --> F["Initialize value = singular_term(canonical_z)"]
    F --> G["Initialize pole_distance = |canonical_z|"]
    G --> H["For each lattice point ω"]
    H --> I["Compute shifted = canonical_z - ω"]
    I --> J{"shifted too close to 0?"}
    J -->|"yes"| Z1["Return PointTooCloseToLatticePoint"]
    J -->|"no"| K["Update pole_distance = min(pole_distance, |shifted|)"]
    K --> L["Add lattice_term(canonical_z, ω) into value"]
    L --> H

    H --> M{"value has finite real and imaginary parts?"}
    M -->|"no"| Z2["Return NumericalComparisonFailed"]
    M -->|"yes"| N["Call build_approximation(z, value, truncation, number_of_terms, pole_distance)"]
    N --> O["Return the caller-chosen approximation/report"]
```
