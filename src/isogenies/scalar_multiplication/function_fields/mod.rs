//! Function-field-side scalar-multiplication story.
//!
//! This submodule keeps two related but distinct narratives separate:
//!
//! - `pullback.rs`: direct construction of `[n]^*` from the generic point
//! - `verschiebung.rs`: the characteristic-`p` factorization
//!   `[p] = V \circ Frob_p` and its certified pullback consequences

mod pullback;
mod verschiebung;
