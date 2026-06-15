//! Structured reports for duality stories between isogenies.
//!
//! This submodule keeps three small narratives separate:
//!
//! - `summaries.rs`: reusable compact value objects
//! - `report.rs`: the main public dual-report surface
//! - `build.rs`: crate-private assembly helpers

mod build;
mod report;
mod summaries;

#[cfg(test)]
mod tests;

pub(crate) use build::DualIsogenySideSummary;
pub use report::DualIsogenyReport;
pub use summaries::{DegreeFactorizationSummary, DualityKind, KernelDescriptionSummary};
