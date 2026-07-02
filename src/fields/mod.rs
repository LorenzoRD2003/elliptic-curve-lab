//! Field-oriented abstractions and tentative data structures.
//!
//! This module is the first implementation target of the project, so the
//! public API is intentionally kept small and documented.

pub mod binary_prime_field;
mod characteristic;
pub mod complex_approx;
pub mod error;
pub mod extension_field;
pub mod finite_field_descriptor;
pub mod montgomery_prime_field;
pub mod polynomial_field;
pub mod rational_function_field;
pub mod rationals;
pub mod traits;

#[cfg(test)]
mod crypto_bigint_spike;

pub use binary_prime_field::{Fp2, Fp2Elem};
pub use characteristic::FieldCharacteristic;
pub use complex_approx::ComplexApprox;
pub use error::FieldError;
pub use montgomery_prime_field::{
    Fp, Fp3, Fp3Elem, Fp5, Fp5Elem, Fp7, Fp7Elem, Fp11, Fp11Elem, Fp13, Fp13Elem, Fp17, Fp17Elem,
    Fp19, Fp19Elem, Fp23, Fp23Elem, Fp29, Fp29Elem, Fp31, Fp31Elem, Fp37, Fp37Elem, Fp41, Fp41Elem,
    Fp43, Fp43Elem, Fp89, Fp89Elem, Fp101, Fp101Elem, Fp241, Fp241Elem, Fp1000000007,
    Fp1000000007Elem, FpElem,
};
pub use rationals::Q;
