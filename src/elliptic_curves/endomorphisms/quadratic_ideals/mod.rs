//! Quadratic-ideal vocabulary for endomorphism-ring computations.
//!
//! This module starts one step before actual ideals: it records how an odd
//! prime `ℓ` behaves in an imaginary quadratic order `O_f`. That local datum is
//! the arithmetic input needed before we can honestly introduce prime ideals,
//! ideal classes, or a class-group action on elliptic curves.
//!
//! The intended staged path is:
//!
//! 1. classify prime behavior in `O_f`;
//! 2. represent small ideals in imaginary quadratic orders;
//! 3. distinguish ideals of prime norm `ℓ`;
//! 4. connect those local prime-norm ideals to horizontal `ℓ`-isogenies;
//! 5. only later promote that local story to a genuine class-group action.
//!
//! In the eventual theory, an ideal class acts on an elliptic curve by
//! `[𝔞] * E = E / E[𝔞]`. This module should own the ideal-side data needed to
//! make that statement precise, while curve quotients and isogeny edges remain
//! owned by the isogeny layers.
//!
//! The current executable surface is intentionally austere: for odd primes not
//! dividing the conductor, [`ImaginaryQuadraticOrder::prime_behavior`] answers
//! the local symbol `(Δ/ℓ)` as split, inert, or ramified data. If `ℓ | f`, the
//! method reports that the prime is not invertible in the non-maximal order
//! instead of pretending that a horizontal ideal action is available.
//!
//! [`ImaginaryQuadraticOrder::prime_behavior`]: super::quadratic_orders::ImaginaryQuadraticOrder::prime_behavior

mod error;
mod prime_behavior;

#[cfg(test)]
mod tests;

pub use error::QuadraticPrimeBehaviorError;
pub use prime_behavior::QuadraticPrimeBehavior;
