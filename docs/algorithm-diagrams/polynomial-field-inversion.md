# Polynomial-Field Inversion

Source: [src/fields/polynomial_field.rs](../../src/fields/polynomial_field.rs)

This follows the same Euclidean idea as the static extension-field path, but
the value carries its own modulus. That makes the distinction between

- a general quotient algebra element, and
- an element of a true quotient field when the modulus is irreducible

more visible at the value layer.

```mermaid
flowchart TB
    A["Input: PolynomialFieldElement with representative f(x) and stored modulus m(x)"] --> B["Reduce f(x) modulo m(x)"]
    B --> C{"Reduced representative is zero?"}
    C -->|"yes"| Z0["Return DivisionByZero"]
    C -->|"no"| D["Run extended Euclidean algorithm on reduced f(x) and m(x)"]
    D --> E["Obtain gcd(x), Bézout coefficient s(x), and t(x)"]
    E --> F{"gcd is a non-zero constant unit?"}
    F -->|"no"| Z1["Return NonInvertibleElement"]
    F -->|"yes"| G["Invert the unit with F::inverse(unit)"]
    G --> H["Scale s(x) by the unit inverse"]
    H --> I["Build a new quotient element from that polynomial"]
    I --> J["Reduce again to the canonical remainder representative"]
    J --> K["Return the inverse in the stored quotient algebra"]

    A --> L["Optional separate check_field_conditions()"]
    L --> M["Ask whether the stored modulus is suitable for a true quotient field"]
    M --> N["If irreducible, every non-zero class should become invertible"]
    M --> O["If reducible, non-zero zero-divisors may still fail at step F"]
```
