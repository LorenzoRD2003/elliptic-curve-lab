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
pub mod traits;
pub mod utils;

pub use crate::visualization::Visualizable;
pub use crate::visualization::fields::{
    VisualizableField, addition_table, describe_complex,
    describe_complex_polynomial_modulus_as_field_modulus, describe_prime_polynomial_field_element,
    describe_prime_polynomial_modulus, describe_prime_polynomial_modulus_as_field_modulus,
    describe_rational, explain_add, explain_complex_polynomial_modulus_irreducibility,
    explain_inverse, explain_mul, explain_prime_polynomial_modulus_irreducibility,
    explain_prime_polynomial_storage, explain_rational_add, explain_rational_div,
    explain_rational_inverse, explain_rational_mul, format_complex, format_complex_polynomial,
    format_fp_elem, format_prime_field, format_prime_polynomial,
    format_prime_polynomial_field_element, format_prime_polynomial_modulus, format_rational,
    format_rational_field, inverses_table, multiplication_table,
};
pub use complex_approx::ComplexApprox;
pub use errors::FieldError;
pub use extension_field::{ExtensionField, ExtensionFieldDescriptor, ExtensionFieldElement};
pub use finite_field::FiniteFieldDescriptor;
pub use polynomial_field::{PolynomialFieldElement, PolynomialModulus};
pub use prime_field::{Fp, FpElem};
pub use rationals::Q;
pub use traits::{Field, FiniteField};
