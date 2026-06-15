//! Irreducibility classification for dense univariate polynomials.
//!
//! This module intentionally owns:
//!
//! - structured result types such as [`IrreducibilityStatus`]
//! - the backend capability [`IrreducibilityBackend`]
//! - backend-specific implementations hidden behind that capability
//!
//! The actual query lives on [`crate::polynomials::DensePolynomial`] through
//! methods such as `polynomial.irreducibility_status()`. This keeps the
//! question attached to the polynomial object while still letting the backend
//! algorithms evolve independently for different coefficient-field families.

mod algebraically_closed;
mod backend;
mod prime_fields;
mod rationals;
mod status;

pub use backend::IrreducibilityBackend;
pub use status::{IrreducibilityStatus, ReducibilityReason};

#[cfg(test)]
mod tests;
