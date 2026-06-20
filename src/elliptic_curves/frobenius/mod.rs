//! The current module keeps two distinctions explicit:
//!
//! - the absolute Frobenius `π_p`, which raises coordinates to `p`-power
//! - the relative Frobenius `π_q`, where `q = p^r` is the size of the
//!   chosen finite base field
//!
//! For a curve defined over `F_q`, `π_q` is an endomorphism of that curve.
//! By contrast, `π_p` need only land on the `p`-power Frobenius twist in
//! general. When the curve coefficients already lie in the prime field, those
//! two curve models coincide, but the fixed points can still differ: this is
//! the first visible distinction between `E(F_p)` and larger finite-field
//! point sets represented in the same coordinate field.
//!
//! In the current educational implementation, both maps act only on point
//! coordinates already represented in one concrete finite field backend. The
//! module does not yet introduce a separate point type for geometric points
//! over an algebraic closure. Even so, the fixed-point story is already
//! visible:
//!
//! - points fixed by `π_p` behave like `E(F_p)` inside the chosen backend
//! - points fixed only by `π_q` can witness larger rationality fields such
//!   as `E(F_{p^r})`

pub mod character_sum;
pub mod characteristic_equation;
pub mod extension_counts;
pub mod group_order;
pub mod hasse;
mod invariants;
mod metadata;
pub mod orbit;
pub mod quadratic_twist;
pub mod schoof;
pub mod torsion;

#[cfg(test)]
mod tests;

pub use crate::elliptic_curves::traits::FrobeniusTraceCurveModel;
pub use hasse::HasseInterval;
pub use invariants::{
    FrobeniusCharacteristicPolynomial, FrobeniusCurveType, FrobeniusDiscriminant,
    FrobeniusLocalZetaFunction, FrobeniusTrace,
};
pub use metadata::{AbsoluteFrobenius, RelativeFrobenius};
