# Lagrange Interpolation

Source: [src/polynomials/interpolation.rs](../../../src/polynomials/interpolation.rs)

This is the current direct educational interpolation routine: for each sample
point, build one Lagrange basis polynomial explicitly, scale it by `y_i`, and
accumulate it into the final result.

```mermaid
flowchart TB
    A["Input: samples (x_0, y_0), ..., (x_{n-1}, y_{n-1})"] --> B{"samples is empty?"}
    B -->|"yes"| Z0["Return the zero polynomial"]
    B -->|"no"| C["Initialize result = 0"]

    C --> D["For each sample index i"]
    D --> E["Initialize numerator = 1 and denominator = 1"]
    E --> F["Loop over every sample index j"]
    F --> G{"j = i?"}
    G -->|"yes"| H["Skip this j"]
    G -->|"no"| I{"x_i = x_j?"}
    I -->|"yes"| Z1["Return DuplicateInterpolationAbscissa"]
    I -->|"no"| J["Multiply numerator by (x - x_j)"]
    J --> K["Multiply denominator by (x_i - x_j)"]
    K --> F
    H --> F

    F --> L["After the inner loop, compute scaling = y_i / denominator"]
    L --> M{"denominator invertible?"}
    M -->|"no"| Z2["Return NonInvertibleInterpolationDenominator"]
    M -->|"yes"| N["Scale the basis polynomial numerator by scaling"]
    N --> O["Add the scaled basis into result"]
    O --> P{"More sample indices i?"}
    P -->|"yes"| D
    P -->|"no"| Q["Return the accumulated polynomial"]
```
