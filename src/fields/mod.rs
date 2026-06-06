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

pub use crate::define_fp_quadratic_extension;
pub use crate::define_q_quadratic_extension;
pub use complex_approx::ComplexApprox;
pub use errors::FieldError;
pub use extension_field::{ExtensionField, ExtensionFieldElement, ExtensionFieldSpec};
pub use finite_field::FiniteFieldDescriptor;
pub use polynomial_field::{PolynomialFieldElement, PolynomialModulus};
pub use prime_field::{Fp, FpElem};
pub use rationals::Q;
pub use sqrt_field::SqrtField;
pub use traits::{EnumerableFiniteField, Field, FiniteField};
