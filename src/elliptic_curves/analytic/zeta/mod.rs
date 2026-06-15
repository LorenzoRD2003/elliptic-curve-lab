//! Truncated Weierstrass `ζ` attached to one lattice.
//!
//! The approximation report stays public in this namespace, while the actual
//! evaluation route is owned by `ComplexLattice`.

mod api;
mod value;

#[cfg(test)]
mod tests;

pub use value::WeierstrassZetaApprox;
