//! Formal composition of explicit isogenies.
//!
//! This module keeps four small stories separate:
//!
//! - `bridge.rs`: how the middle curves are matched or transported
//! - `composed_isogeny.rs`: the public composed-isogeny type and constructors
//! - `evaluation.rs`: pointwise evaluation of `second ∘ bridge ∘ first`
//! - `kernel.rs`: exhaustive kernel enumeration for the composed map

mod bridge;
mod composed_isogeny;
mod evaluation;
mod kernel;

#[cfg(test)]
mod tests;

pub use composed_isogeny::ComposedIsogeny;
