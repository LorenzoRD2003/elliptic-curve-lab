//! Foundational scaffolding for mathematical and cryptographic algorithms.
//!
//! The crate intentionally starts with small, documented interfaces and
//! lightweight placeholder implementations so the core abstractions can evolve
//! with tests before the heavy algebraic algorithms arrive.

// pub mod algorithms;
// pub mod elliptic_curves;
pub mod fields;
pub mod polynomials;
pub mod visualization;

// pub use elliptic_curves::{AffinePoint, CurveEquation, ProjectivePoint, ShortWeierstrassCurve};
pub use fields::{
    ComplexApprox, ExtensionField, ExtensionFieldElement, Field, FieldError, FiniteField,
    FiniteFieldDescriptor, Fp, FpElem, PolynomialFieldElement, PolynomialModulus, Q,
    addition_table, describe_complex, describe_complex_polynomial_modulus_as_field_modulus,
    describe_prime_polynomial_field_element, describe_prime_polynomial_modulus,
    describe_prime_polynomial_modulus_as_field_modulus, describe_rational, explain_add,
    explain_complex_polynomial_modulus_irreducibility, explain_inverse, explain_mul,
    explain_prime_polynomial_modulus_irreducibility, explain_prime_polynomial_storage,
    explain_rational_add, explain_rational_div, explain_rational_inverse, explain_rational_mul,
    format_complex, format_complex_polynomial, format_fp_elem, format_prime_field,
    format_prime_polynomial, format_prime_polynomial_field_element,
    format_prime_polynomial_modulus, format_rational, format_rational_field, inverses_table,
    multiplication_table,
};
pub use polynomials::{
    DensePolynomial, IrreducibilityBackend, IrreducibilityStatus, PolynomialError,
    ReducibilityReason, SparsePolynomial, VisualizablePolynomial, describe_dense_polynomial,
    describe_multivariate_polynomial, describe_sparse_polynomial, explain_dense_storage,
    explain_multivariate_storage, explain_sparse_storage, format_dense_polynomial, format_monomial,
    format_multivariate_polynomial, format_sparse_polynomial, irreducibility_status,
    is_irreducible,
};
pub use visualization::{Visualizable, VisualizableField};
