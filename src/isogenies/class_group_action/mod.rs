//! Introductory scaffolding for class-group actions by horizontal isogenies.
//!
//! This module names the bridge between ideal classes in imaginary quadratic
//! orders and isogeny-graph motion. The library already implements finite
//! quadratic class groups algebraically through reduced binary quadratic forms,
//! composition, subgroups, and Cayley tables; this module does not yet
//! implement the full arithmetic CM action through ideal kernels.
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

mod action_plan;
mod crater_walk;
mod graph_reports;
mod horizontal_ideal;
mod ideal_label;
mod isogeny_action;
#[allow(dead_code)]
mod labeled_crater_walk;
mod orientation;
mod witness_set;

#[cfg(test)]
mod tests;

pub use action_plan::{ClassGroupActionPlan, ClassGroupActionPlanError};
pub use crater_walk::{CraterWalkReport, CraterWalkTermination};
pub use horizontal_ideal::{HorizontalIdealReport, HorizontalIdealStatus};
pub use ideal_label::CraterIdealLabelError;
pub use isogeny_action::{
    ClassGroupIsogenyActionError, ClassGroupIsogenyActionReport, ClassGroupIsogenyActionSegment,
};
pub use labeled_crater_walk::{
    CraterDirectionCertification, LabeledCraterWalkError, LabeledCraterWalkReport,
};
pub use orientation::{
    CraterOrientationWitness, CraterOrientationWitnessError, OrientedCraterClassOrderComparison,
    OrientedCraterClassOrderComparisonError, OrientedCraterClassOrderStatus,
    OrientedCraterPowerActionError, OrientedCraterPowerActionReport,
    OrientedLabeledCraterWalkReport,
};

pub(crate) use action_plan::ClassGroupActionPlanFactor;
pub(crate) use horizontal_ideal::HorizontalIdealWitness;
pub(crate) use ideal_label::{CraterIdealLabelReport, CraterIdealPrimeBehavior};
pub(crate) use witness_set::LocalCraterWitnessSet;
