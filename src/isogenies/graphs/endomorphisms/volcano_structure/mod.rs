//! Certified structural views of a stored ordinary `ℓ`-volcano graph.
//!
//! This module builds a global report from the same floor-distance evidence as
//! the existing shortest-floor search. It treats
//! `δ(v) = dist(v, V_d)` as an altimeter and reconstructs certified levels by
//! subtracting each distance from the largest certified distance seen in the
//! stored graph.
//!
//! The resulting structure is intentionally an evidence report for the stored
//! graph. If some vertices are partial, special, or otherwise fail the
//! ordinary-volcano floor search, they remain outside the certified level
//! partition and are reported separately.

mod build;
mod level;
mod node;
mod report;
mod uncertified;

#[cfg(test)]
mod tests;

pub use level::{VolcanoStructureLevelReport, VolcanoStructureRole};
pub use node::VolcanoStructureNodeReport;
pub use report::VolcanoStructureReport;
pub use uncertified::{UncertifiedVolcanoNodeReport, VolcanoStructureUncertifiedReason};
