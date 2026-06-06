# Extension-Field Inversion

Source: [src/fields/extension_field.rs](../../../src/fields/extension_field.rs)

This is the type-level field-family version: the ambient modulus lives in the
`ExtensionFieldSpec`, while each element stores only its quotient
representative.

```mermaid
flowchart TB
    A["Input: extension-field element [f(x)] in Base[x]/(m(x))"] --> B["Reduce f(x) modulo m(x)"]
    B --> C{"Reduced representative is zero?"}
    C -->|"yes"| Z0["Return DivisionByZero"]
    C -->|"no"| D["Set representative = reduced f(x), modulus = m(x)"]
    D --> E["Run extended Euclidean algorithm on representative and modulus"]
    E --> F["Obtain gcd(x), Bézout coefficient s(x), and t(x)"]
    F --> G{"gcd is a non-zero constant unit?"}
    G -->|"no: no constant term or degree > 0"| Z1["Return NonInvertibleElement"]
    G -->|"yes"| H["Invert the unit in the base field"]
    H --> I["Scale the Bézout coefficient by unit_inverse"]
    I --> J["Reduce the scaled Bézout polynomial modulo m(x)"]
    J --> K["Return the canonical inverse class"]
```
