//! Floor-search evidence for ordinary `ℓ`-isogeny volcanoes.
//!
//! This module implements the first, deliberately modest, part of §3.1 of
//! Sutherland's "Isogeny volcanoes": use local outgoing degree evidence to
//! recognize floor vertices in an ordinary `ℓ`-volcano. The API stays honest
//! about the current graph builder: a stored node may be present without having
//! its outgoing `ℓ`-isogenies fully expanded, so low observed degree is only
//! decisive when the node is known to be complete.

mod error;
mod evidence;
mod path;

#[cfg(test)]
mod tests;

pub use error::VolcanoSearchError;
pub use path::FloorPathReport;
