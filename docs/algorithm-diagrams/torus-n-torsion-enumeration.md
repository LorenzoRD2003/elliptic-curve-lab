# Torus N-Torsion Enumeration

Source: [src/elliptic_curves/analytic/torsion/torus.rs](../../src/elliptic_curves/analytic/torsion/torus.rs)

These helpers materialize the reduced `n × n` torsion grid in lexicographic
order, and optionally filter it to the primitive classes of exact torus order
`n`.

```mermaid
flowchart TB
    A["Input: lattice Λ and order n"] --> B{"n = 0?"}
    B -->|"yes"| Z0["Return InvalidTorusTorsionIndex"]
    B -->|"no"| C["Initialize points = []"]
    C --> D["Loop over a = 0 .. n-1"]
    D --> E["Loop over b = 0 .. n-1"]
    E --> F["Build reduced torsion index (a, b; n)"]
    F --> G["Convert it to the fundamental coordinate (a/n, b/n)"]
    G --> H["Map that coordinate to z = (a/n)ω₁ + (b/n)ω₂"]
    H --> I["Store TorusTorsionPoint { index, coordinate, z }"]
    I --> E
    E --> J{"More a values?"}
    J -->|"yes"| D
    J -->|"no"| K["Return the full torus n-torsion list"]

    K --> L["primitive_torus_n_torsion_points"]
    L --> M["Filter the full list by index.is_primitive()"]
    M --> N["Return only the primitive torus n-torsion classes"]
```
