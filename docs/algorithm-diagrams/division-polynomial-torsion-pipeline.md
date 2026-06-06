# Division-Polynomial Torsion Pipeline

Source: [src/elliptic_curves/division_polynomials/torsion.rs](../../src/elliptic_curves/division_polynomials/torsion.rs)

This is the main pedagogical division-polynomial pipeline: start from an index `n`,
dispatch by parity, build rational candidates, then refine them into actual
`n`-torsion and exact-order-`n` points, with an optional comparison against
exhaustive enumeration.

```mermaid
flowchart TB
    A["Input: curve E and index n"] --> B{"n = 0?"}
    B -->|"yes"| Z0["Return ZeroIndex error"]
    B -->|"no"| C["Classify parity with division_polynomial_x_criterion_kind(n)"]

    C -->|"odd n"| D["Build odd division polynomial ψ_n(x)"]
    D --> E["Enumerate x in F and keep roots of ψ_n(x)"]
    E --> F["Keep only x that lift to at least one rational affine point"]
    F --> G["Lift each surviving x to affine candidates P = (x, y)"]

    C -->|"even n"| H["Use the stripped even factor ε_n(x) in ψ_n = y ε_n(x)"]
    H --> I["Enumerate rational affine points P = (x, y) on E"]
    I --> J{"y(P) = 0 or ε_n(x(P)) = 0?"}
    J -->|"yes"| K["Keep P as a raw even-case candidate"]
    J -->|"no"| L["Discard P"]

    G --> M["torsion_candidates_from_division_polynomial"]
    K --> M

    M --> N{"n odd?"}
    N -->|"yes"| O["Every candidate is kept as n-torsion"]
    N -->|"no"| P["Retain y = 0 candidates only when [n]P = O"]

    O --> Q["torsion_points_from_division_polynomial"]
    P --> Q

    Q --> R["Check exact order point-by-point with point_has_exact_order"]
    R --> S["exact_n_torsion_points_from_division_polynomial"]

    S --> T["Optionally compare against exhaustive enumeration of E(F_q)"]
    T --> U["Report counts, missing exact-order points, and extras"]
```
