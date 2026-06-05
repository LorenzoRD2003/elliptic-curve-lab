# Isogeny Graph Scaffolding

Source: [src/isogenies/graphs](../../src/isogenies/graphs)

This note was previously stored as `src/isogenies/graphs/README.md`.

This directory contains the educational `\ell`-isogeny graph scaffolding for
small short-Weierstrass curves over prime fields.

The key modeling choice is that each node stores one chosen representative
curve, while each edge stores:

- a rational cyclic kernel on the source representative
- the directed source and target node ids
- an optional witness transporting the raw Vélu codomain onto the stored
  target representative

That separation keeps two different notions visible:

- the codomain curve produced directly by Vélu from a kernel
- the representative curve chosen for the target node after deduplication

```mermaid
flowchart LR
    A["Source node<br/>representative E"] --> B["Pick rational cyclic kernel G = &lt;P&gt;"]
    B --> C["Build Vélu isogeny φ_raw : E → E/G"]
    C --> D{"Does E/G already match<br/>a stored representative?"}
    D -->|"yes, exact match"| E["Use existing target node<br/>witness = Identity"]
    D -->|"yes, up to base-field isomorphism"| F["Use existing target node<br/>witness = explicit isomorphism α"]
    D -->|"no"| G["Insert fresh target node<br/>representative = E/G"]

    E --> H["Store graph edge<br/>α ∘ φ_raw collapses to φ_raw"]
    F --> I["Store graph edge<br/>map interpreted as α ∘ φ_raw"]
    G --> J["Store graph edge<br/>witness = Identity"]

    H --> K["Local verification"]
    I --> K
    J --> K

    K --> L["Check map lands on target representative"]
    K --> M["Check kernel maps to identity"]
    K --> N["Check homomorphism law on rational points"]
    K --> O["Inspect reverse-direction edges already in graph"]
    O --> P["Classify as Missing / PresentButNotVerifiedAsDual / VerifiedAsDual"]
```

In other words, the graph stores representatives and witnesses explicitly so
that later summaries and local verification can reason about the maps that
actually connect the stored nodes, not just about abstract isomorphism classes.
