//! Core trait surfaces for explicit isogenies.
//!
//! This submodule separates:
//!
//! - `isogeny.rs`: the minimal shared map/kernel interface
//! - `degree_factorized.rs`: separable/inseparable degree metadata
//! - `verifiable.rs`: exhaustive small-field verification helpers

mod degree_factorized;
mod isogeny;
mod verifiable;

#[cfg(test)]
mod tests;

pub use degree_factorized::DegreeFactorizedIsogeny;
pub use isogeny::Isogeny;
pub use verifiable::VerifiableIsogeny;
