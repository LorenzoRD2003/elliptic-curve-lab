# Prime-Field Square Root

Source: [src/fields/prime_field.rs](../../src/fields/prime_field.rs)

This is the current exact prime-field square-root routine. It handles the
small easy cases directly, rejects non-residues honestly, uses the
`p % 4 == 3` shortcut when available, and otherwise falls back to the full
Tonelli–Shanks loop.

```mermaid
flowchart TB
    A["Input: x in Fp(P)"] --> B["Validate the prime-field modulus P"]
    B --> C{"modulus valid?"}
    C -->|"no"| Z0["Return None"]
    C -->|"yes"| D{"x = 0?"}
    D -->|"yes"| E["Return 0"]
    D -->|"no"| F{"P = 2?"}
    F -->|"yes"| G["Return x itself"]
    F -->|"no"| H{"x is a quadratic residue?"}
    H -->|"no"| Z1["Return None"]
    H -->|"yes"| I{"P % 4 = 3?"}
    I -->|"yes"| J["Return x^((P + 1) / 4)"]
    I -->|"no"| K["Decompose P - 1 = q * 2^s with q odd"]
    K --> L["Find a quadratic non-residue z"]
    L --> M["Initialize m = s, c = z^q, t = x^q, r = x^ceil(q/2)"]
    M --> N{"t = 1?"}
    N -->|"yes"| O["Return r"]
    N -->|"no"| P["Find the least i < m with t^(2^i) = 1"]
    P --> Q{"such i exists?"}
    Q -->|"no"| Z2["Return None"]
    Q -->|"yes"| R["Set b = c^(2^(m - i - 1))"]
    R --> S["Update r = r*b, t = t*b^2, c = b^2, m = i"]
    S --> N
```
