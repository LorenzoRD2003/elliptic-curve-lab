//! Field-oriented abstractions and tentative data structures.
//!
//! This module is the first implementation target of the project, so the
//! public API is intentionally kept small and documented.

pub mod complex_approx;
pub mod error;
pub mod extension_field;
pub mod finite_field_descriptor;
pub mod polynomial_field;
pub mod prime_field;
pub mod rational_function_field;
pub mod rationals;
pub mod traits;

#[cfg(test)]
mod crypto_bigint_spike;

pub use complex_approx::ComplexApprox;
pub use error::FieldError;
pub use prime_field::{Fp, FpElem};
pub use rationals::Q;
