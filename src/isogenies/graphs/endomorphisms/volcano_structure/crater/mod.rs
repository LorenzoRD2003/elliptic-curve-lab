//! Crater and horizontal-edge evidence for a stored `ℓ`-volcano graph.
//!
//! The crater is the certified surface `V_0` of a
//! [`VolcanoStructureReport`](super::VolcanoStructureReport). This submodule
//! records the horizontal edges observed on that surface, while keeping
//! certified altimeter evidence separate from weaker graph-theoretic hints.

mod build;
mod edge;
mod report;
mod shape;

#[cfg(test)]
mod tests;

pub use edge::{HorizontalEdgeReport, HorizontalEdgeStatus};
pub use report::CraterReport;
pub use shape::CraterShape;
