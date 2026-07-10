//! Introductory scaffolding for class-group actions by horizontal isogenies.
//!
//! This module names the bridge between ideal classes in imaginary quadratic
//! orders and isogeny-graph motion, without implementing a full class group or
//! a general `E[a]` kernel construction yet.
//!
//! The guiding mathematical slogan is:
//!
//! ```text
//! [a] * E = E / E[a]
//! ```
//!
//! The first executable milestone should be much smaller than that general
//! formula: a local prime-norm story where an ideal of norm `ell` corresponds
//! to a horizontal edge in an `ell`-isogeny volcano. That keeps the current
//! volcano evidence graph-structural until the arithmetic layer can certify the
//! ideal-side interpretation.
//!
//! Ownership is split deliberately:
//!
//! - `elliptic_curves::endomorphisms::quadratic_ideals` should own ideal-side
//!   objects and local ideal-norm vocabulary;
//! - this module should own reports and adapters that interpret certified
//!   horizontal `ell`-isogeny evidence as the first shadow of an ideal action.

mod crater_walk;
mod graph_reports;
mod horizontal_ideal;
mod ideal_label;
#[allow(dead_code)]
mod labeled_crater_walk;
mod orientation;

#[cfg(test)]
mod tests;

pub use crater_walk::{CraterWalkReport, CraterWalkTermination};
pub use horizontal_ideal::{HorizontalIdealReport, HorizontalIdealStatus, HorizontalIdealWitness};
pub use ideal_label::{CraterIdealLabelError, CraterIdealLabelReport, CraterIdealPrimeBehavior};
pub use labeled_crater_walk::{
    CraterDirectionCertification, LabeledCraterWalkError, LabeledCraterWalkReport,
};
pub use orientation::{
    CraterOrientationWitness, CraterOrientationWitnessError, OrientedCraterClassOrderComparison,
    OrientedCraterClassOrderComparisonError, OrientedCraterClassOrderStatus,
    OrientedCraterPowerActionError, OrientedCraterPowerActionReport,
    OrientedLabeledCraterWalkReport,
};
