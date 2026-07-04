//! Imaginary quadratic-order arithmetic for the current `endomorphisms` layer.
//!
//! The implemented public story in this subtree is: integral quadratic
//! discriminants, their canonical factorization `Δ = f^2 D_K`, and
//! imaginary quadratic orders `O_f = ℤ + f O_K`.
//!
//! We want to keep the pipeline aligned with the consumers: `m -> D_K -> O_K`.

mod cover_relation;
mod discriminant;
mod error;
mod factorization;
mod order;

pub use cover_relation::QuadraticOrderCoverRelation;
pub use discriminant::{DiscriminantSign, QuadraticDiscriminant};
pub use error::{
    ImaginaryQuadraticOrderError, QuadraticDiscriminantFactorizationError, QuadraticOrderIndexError,
    QuadraticRadicandError,
};
pub use factorization::QuadraticDiscriminantFactorization;
pub use order::ImaginaryQuadraticOrder;
