//! Foundational scaffolding for mathematical and cryptographic algorithms.
//!
//! The crate intentionally starts with small, documented interfaces and
//! lightweight placeholder implementations so the core abstractions can evolve
//! with tests before the heavy algebraic algorithms arrive.
//!
//! The crate root is intentionally just a module index.
//! Public entry points live under their mathematical namespaces such as
//! `fields`, `polynomials`, `elliptic_curves`, `isogenies`, and `numerics`.

pub mod elliptic_curves;
pub mod fields;
pub mod isogenies;
pub mod numerics;
pub mod polynomials;
#[cfg(feature = "visualization")]
pub mod visualization;

// Internal property-testing scaffolding. This stays crate-private even when
// the `test-support` feature is enabled; the feature only widens internal
// availability, not the public API surface.
#[cfg(any(test, feature = "test-support"))]
pub(crate) mod proptest_support;
