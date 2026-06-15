mod cover_relation;
mod discriminant;
mod error;
mod factorization;
mod order;

pub use cover_relation::QuadraticOrderCoverRelation;
pub use discriminant::{DiscriminantSign, QuadraticDiscriminant};
pub use error::{
    ImaginaryQuadraticOrderError, QuadraticDiscriminantFactorizationError, QuadraticOrderIndexError,
};
pub use factorization::QuadraticDiscriminantFactorization;
pub use order::ImaginaryQuadraticOrder;

pub(crate) use discriminant::QuadraticDiscriminantMod4;
