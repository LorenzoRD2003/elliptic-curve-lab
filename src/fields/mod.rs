//! Field-oriented abstractions and tentative data structures.
//!
//! This module is the first implementation target of the project, so the
//! public API is intentionally kept small and documented.

pub mod complex_approx;
pub mod errors;
pub mod extension_field;
pub mod finite_field;
pub mod polynomial_field;
pub mod prime_field;
pub mod rationals;
pub mod sqrt_field;
pub mod traits;
pub mod utils;

pub use crate::visualization::Visualizable;
pub use crate::visualization::fields::{
    VisualizableField, addition_table, describe_complex,
    describe_complex_polynomial_modulus_as_field_modulus, describe_extension_field,
    describe_extension_field_element, describe_prime_polynomial_field_element,
    describe_prime_polynomial_modulus, describe_prime_polynomial_modulus_as_field_modulus,
    describe_rational, explain_add, explain_complex_polynomial_modulus_irreducibility,
    explain_complex_square_root, explain_extension_field_add, explain_extension_field_inverse,
    explain_extension_field_mul, explain_extension_field_reduction, explain_inverse, explain_mul,
    explain_prime_field_square_root, explain_prime_polynomial_field_add,
    explain_prime_polynomial_field_inverse, explain_prime_polynomial_field_mul,
    explain_prime_polynomial_field_reduction, explain_prime_polynomial_modulus_irreducibility,
    explain_prime_polynomial_storage, explain_rational_add, explain_rational_div,
    explain_rational_inverse, explain_rational_mul, explain_rational_square_root, format_complex,
    format_complex_polynomial, format_extension_field, format_extension_field_element,
    format_fp_elem, format_prime_field, format_prime_polynomial,
    format_prime_polynomial_field_element, format_prime_polynomial_modulus, format_rational,
    format_rational_field, inverses_table, multiplication_table,
};
pub use complex_approx::ComplexApprox;
pub use errors::FieldError;
pub use extension_field::{ExtensionField, ExtensionFieldElement, ExtensionFieldSpec};
pub use finite_field::FiniteFieldDescriptor;
pub use polynomial_field::{PolynomialFieldElement, PolynomialModulus};
pub use prime_field::{Fp, FpElem};
pub use rationals::Q;
pub use sqrt_field::SqrtField;
pub use traits::{EnumerableFiniteField, Field, FiniteField};
