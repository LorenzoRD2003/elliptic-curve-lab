//! Foundational scaffolding for mathematical and cryptographic algorithms.
//!
//! The crate intentionally starts with small, documented interfaces and
//! lightweight placeholder implementations so the core abstractions can evolve
//! with tests before the heavy algebraic algorithms arrive.

// pub mod algorithms;
// pub mod elliptic_curves;
pub mod fields;
// pub mod polynomials;

// pub use elliptic_curves::{AffinePoint, CurveEquation, ProjectivePoint, ShortWeierstrassCurve};
pub use fields::{
    ComplexApprox, ExtensionField, ExtensionFieldElement, Field, FieldError, FiniteField,
    FiniteFieldDescriptor, Fp, FpElem, PolynomialFieldElement, PolynomialModulus, Q, Visualizable,
    addition_table, describe_complex, describe_prime_polynomial_field_element,
    describe_prime_polynomial_modulus, describe_rational, explain_add, explain_inverse,
    explain_mul, explain_prime_polynomial_storage, explain_rational_add, explain_rational_div,
    explain_rational_inverse, explain_rational_mul, format_complex, format_fp_elem,
    format_prime_field, format_prime_polynomial, format_prime_polynomial_field_element,
    format_prime_polynomial_modulus, format_rational, format_rational_field, inverses_table,
    multiplication_table,
};
// pub use polynomials::{DensePolynomial, SparsePolynomial};
