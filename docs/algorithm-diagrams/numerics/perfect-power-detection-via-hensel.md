# Perfect-Power Detection Via Hensel

Source:
[src/numerics/perfect_powers/mod.rs](../../../src/numerics/perfect_powers/mod.rs)

This staged detector handles positive integers `N` with `gcd(N, 6) = 1`. It
tests prime exponents `q ≤ ⌊log₂ N⌋` by converting the question `N = aᵠ` into
the integer-root problem

`f_q(x) = xᵠ − N`.

For each exponent, it uses the integer-root Hensel surface to certify candidate
bases, then performs the final exact check `aᵠ = N`.

```mermaid
flowchart TB
    A["Input: N"] --> B{"N > 1 and gcd(N, 6) = 1?"}
    B -->|"no"| Z0["Return the matching limitation outcome"]
    B -->|"yes"| C["List prime exponents q ≤ ⌊log₂ N⌋"]

    C --> D{"More q?"}
    D -->|"no"| Z1["Return NotPerfectPower"]
    D -->|"yes"| E["Build f_q(x) = xᵠ − N"]

    E --> F{"q = 2?"}
    F -->|"yes"| G["Use Hensel prime p = 3"]
    F -->|"no"| H["Use Hensel prime p = 2"]

    G --> I["Find bounded integer roots of f_q"]
    H --> I
    I --> J["For each certified root r, set a = |r| when positive"]
    J --> K{"a^q = N exactly?"}
    K -->|"yes"| Z2["Return PerfectPower { a, q }"]
    K -->|"no"| D
```

The small-prime choice keeps the simple-root Hensel route honest in this staged
setting. When `N` is coprime to `6`, any genuine base `a` is a unit modulo `2`
and `3`; choosing `p = 3` for `q = 2` and `p = 2` for odd prime `q` ensures
`p ∤ q`, so `q·a^(q−1)` is non-zero modulo `p`.

Complexity: let `n = ⌈log₂ N⌉` and let `M(n)` be the cost of multiplying
`n`-bit integers. The detector tests `π(n) = Θ(n/log n)` prime exponents. For one
exponent `q`, the sparse polynomial `xᵠ − N` has `O(1)` terms, evaluation costs
`O(log q · M(n))`, and the Hensel precision is `Θ(n/q)`. Summing over prime
`q ≤ n` gives `O(n log n · M(n))` bit operations, plus lower-order sieve and
exact-power-check work. With quasi-linear integer multiplication, this is
quasi-quadratic in `n`.
