# Analytic Torsion vs Division Polynomial

Source: [src/elliptic_curves/analytic/torsion/division_polynomial.rs](../../src/elliptic_curves/analytic/torsion/division_polynomial.rs)

This comparison pipeline starts from torus `n`-torsion on `ℂ/Λ`, maps each
class to the analytic curve, and then classifies the division-polynomial
comparison by pole / odd / even case.

```mermaid
flowchart TB
    A["Input: lattice Λ, order n, truncations, tolerance"] --> B["Classify parity with division_polynomial_x_criterion_kind(n)"]
    B --> C["Build analytic Weierstrass curve from Λ"]
    C --> D["Convert it to the short-Weierstrass companion"]
    D --> E["Map every torus n-torsion class to the curve via z ↦ (℘(z), ℘′(z))"]
    E --> F["For each mapped torsion point"]

    F --> G{"Curve point is infinity?"}
    G -->|"yes"| H["Return Pole case with status = PoleAtIdentity"]
    G -->|"no"| I["Read affine coordinates (x, y)"]
    I --> J["Evaluate the complex x-criterion at x"]
    J --> K["Compute absolute value and approximate-zero verdict"]

    K --> L{"Odd criterion ψ_n(x) or even criterion ε_n(x)?"}
    L -->|"odd"| M["Return Odd report"]
    M --> N{"ψ_n(x) approximately zero?"}
    N -->|"yes"| O["Status = VanishesApproximately"]
    N -->|"no"| P["Status = DoesNotVanishApproximately"]

    L -->|"even"| Q["Test whether y is approximately zero"]
    Q --> R["Combine y-branch and ε_n(x)-branch"]
    R --> S["Classify branch as Both / YApproxZero / XCriterionApproxZero / Neither"]
    S --> T{"y approximately zero or ε_n(x) approximately zero?"}
    T -->|"yes"| U["Status = VanishesApproximately"]
    T -->|"no"| V["Status = DoesNotVanishApproximately"]
    U --> W["Return Even report"]
    V --> W

    H --> X["Collect all comparison cases"]
    O --> X
    P --> X
    W --> X
```
