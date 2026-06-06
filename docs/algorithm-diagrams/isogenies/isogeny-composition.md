# Isogeny Composition

Source: [src/isogenies/composition.rs](../../../src/isogenies/composition.rs)

The composition surface supports two pedagogical modes:

- strict composition, where `codomain(first) = domain(second)` exactly
- bridged composition, where an explicit middle isomorphism transports the
  raw middle image onto the second map's chosen domain

```mermaid
flowchart TB
    A["Choose first : E → E' and second : E' or E'' entry map"] --> B{"Strict or bridged?"}

    B -->|"strict"| C["Validate codomain(first) == domain(second)"]
    B -->|"bridged"| D["Provide bridge α : codomain(first) → domain(second)"]
    D --> E["Validate α.domain = codomain(first) and α.codomain = domain(second)"]

    C --> F["Compute rational kernel of the composed map"]
    E --> F

    F --> G["Enumerate every rational point P in domain(first)"]
    G --> H["Evaluate second ∘ bridge ∘ first on P"]
    H --> I{"Image is codomain identity?"}
    I -->|"yes"| J["Include P in kernel_points"]
    I -->|"no"| K["Skip P"]

    J --> L["Construct ComposedIsogeny"]
    K --> L

    L --> M["evaluate(point)"]
    M --> N{"point lies on domain(first)?"}
    N -->|"no"| Z1["Return PointNotOnCurve"]
    N -->|"yes"| O["Compute middle = first(point)"]
    O --> P{"middle lies on codomain(first)?"}
    P -->|"no"| Z2["Return ImagePointNotOnCodomain"]
    P -->|"yes"| Q["Transport middle with identity bridge or α"]
    Q --> R{"bridged point lies on domain(second)?"}
    R -->|"no"| Z3["Return ImagePointNotOnCodomain"]
    R -->|"yes"| S["Compute image = second(bridged point)"]
    S --> T{"image lies on codomain(second)?"}
    T -->|"no"| Z4["Return ImagePointNotOnCodomain"]
    T -->|"yes"| U["Return composed image"]
```
