//! Binary quadratic forms and representation helpers.
//!
//! This module starts with the diagonal positive form `x² + d y²`, which is
//! the representation problem needed by the first CM-oriented workflows.
//! Cornacchia's algorithm is the current engine for primitive representations,
//! but the public surface is phrased in terms of forms rather than in terms of
//! one algorithm.

mod diagonal;
mod error;
mod representation;

#[cfg(test)]
mod tests;

pub use diagonal::DiagonalBinaryQuadraticForm;
pub use error::QuadraticFormError;
pub use representation::DiagonalBinaryQuadraticRepresentation;
