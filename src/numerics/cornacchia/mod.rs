//! Cornacchia's algorithm for equations `x² + d y² = m`.

mod api;
mod error;
mod root;
mod solution;
mod validation;

#[cfg(test)]
mod tests;

pub(crate) use api::cornacchia_candidate_solutions;
pub use api::{cornacchia_primitive_solutions, cornacchia_with_root};
pub use error::CornacchiaError;
pub use solution::CornacchiaSolution;
