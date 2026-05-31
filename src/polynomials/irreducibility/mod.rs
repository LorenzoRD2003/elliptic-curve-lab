//! Irreducibility classification for dense univariate polynomials.
//!
//! The public API of this module is intentionally small:
//!
//! - [`irreducibility_status`] for the structured result
//! - [`is_irreducible`] for the boolean convenience wrapper
//! - backend-specific implementations hidden behind [`IrreducibilityBackend`]
//!
//! This keeps the user-facing surface stable while the internal algorithms can
//! evolve independently for different base-field families.

mod algebraically_closed;
mod prime_fields;
mod rationals;
mod status;
mod traits;

pub use status::{IrreducibilityStatus, ReducibilityReason};
pub use traits::IrreducibilityBackend;

use crate::polynomials::{DensePolynomial, PolynomialError};

/// Returns a structured irreducibility classification for a dense polynomial.
///
/// The exact strategy depends on the base-field backend:
///
/// - prime fields `Fp<P>` currently use an exhaustive educational search over
///   monic candidate divisors
/// - algebraically closed backends such as `ComplexApprox` can conclude that
///   every degree-`>= 2` polynomial is reducible, even when no explicit
///   factorization witness is currently returned
/// - `Q` currently uses an exact but partial backend that either certifies an
///   answer or returns a typed inconclusive error
///
/// Classification conventions:
///
/// - degree `0` polynomials are reported as [`IrreducibilityStatus::Constant`]
/// - degree `1` polynomials are reported as [`IrreducibilityStatus::Linear`]
/// - constants are not considered irreducible
/// - linears are considered irreducible
///
/// TODO:
/// - add a Rabin-style irreducibility test for finite fields once the
///   supporting polynomial infrastructure is mature enough
/// - extend the exact partial backend for `Q` into a complete decision
///   procedure once integer-polynomial factorization infrastructure exists
pub fn irreducibility_status<F: IrreducibilityBackend>(
    polynomial: &DensePolynomial<F>,
) -> Result<IrreducibilityStatus<F>, PolynomialError> {
    F::irreducibility_status_impl(polynomial)
}

/// Returns whether a dense polynomial is irreducible in the current backend.
///
/// This is the boolean convenience wrapper around [`irreducibility_status`].
pub fn is_irreducible<F: IrreducibilityBackend>(
    polynomial: &DensePolynomial<F>,
) -> Result<bool, PolynomialError> {
    Ok(matches!(
        irreducibility_status(polynomial)?,
        IrreducibilityStatus::Linear | IrreducibilityStatus::Irreducible
    ))
}

#[cfg(test)]
mod tests;
