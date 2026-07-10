use core::fmt;

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::QuadraticClassGroup,
    quadratic_ideals::{IdealFormCorrespondence, PrimeNormIdeal},
};
use crate::isogenies::{
    class_group_action::{
        CraterIdealLabelError, CraterIdealLabelReport, CraterOrientationWitness,
        CraterOrientationWitnessError, CraterWalkReport, OrientedLabeledCraterWalkReport,
    },
    graphs::{
        IsogenyGraphNodeId,
        endomorphisms::{CraterReport, VolcanoSearchError},
    },
};

/// Certification status for the direction used by a labeled crater walk.
///
/// The first labeled-walk milestone deliberately records only that the walk
/// direction was chosen deterministically from graph data. It does not certify
/// which horizontal direction corresponds to `𝔭` rather than `\bar{𝔭}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CraterDirectionCertification {
    /// The walk follows the deterministic graph-local rule from [`CraterWalkReport`].
    GraphDeterministic,
    /// A user supplied and graph-validated a crater orientation witness.
    UserSuppliedArithmeticOrientation,
    /// A future implementation inferred the arithmetic crater orientation.
    CertifiedArithmeticOrientation,
}

/// A crater walk equipped with the form class attached to its prime-norm ideal.
///
/// This report connects three independently certified pieces:
///
/// - local ideal/crater/class-group compatibility;
/// - the deterministic graph walk on certified horizontal crater edges;
/// - the reduced binary-quadratic-form label associated to the ideal.
///
/// It still does not claim that the chosen graph direction is the arithmetic
/// action of the form class on the starting curve.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LabeledCraterWalkReport {
    local_label: CraterIdealLabelReport,
    walk: CraterWalkReport,
    form_label: IdealFormCorrespondence,
    direction_certification: CraterDirectionCertification,
}

impl LabeledCraterWalkReport {
    /// Builds a labeled deterministic crater walk from one certified crater report.
    ///
    /// The constructor first checks local compatibility of the supplied ideal
    /// with the crater prime and the class-group discriminant. It then derives
    /// the ideal's binary-quadratic-form class label and runs the existing
    /// deterministic crater walk.
    ///
    /// Complexity: one local compatibility check, one ideal-to-form
    /// conversion, and one deterministic crater walk through certified crater
    /// edges.
    pub(crate) fn from_crater_report(
        crater: &CraterReport,
        class_group: &QuadraticClassGroup,
        ideal: PrimeNormIdeal,
        start: IsogenyGraphNodeId,
    ) -> Result<Self, CraterIdealLabelError> {
        let local_label = CraterIdealLabelReport::new(crater, class_group, ideal.clone())?;
        let form_label = IdealFormCorrespondence::from_prime_norm_ideal(&ideal).expect(
            "a validated prime-norm ideal in the same class group should produce a form label",
        );
        let walk = CraterWalkReport::from_crater_report(crater, ideal, start);

        Ok(Self {
            local_label,
            walk,
            form_label,
            direction_certification: CraterDirectionCertification::GraphDeterministic,
        })
    }

    /// Returns the local ideal/crater/class-group compatibility certificate.
    pub fn local_label(&self) -> &CraterIdealLabelReport {
        &self.local_label
    }

    /// Returns the deterministic crater walk.
    pub fn walk(&self) -> &CraterWalkReport {
        &self.walk
    }

    /// Returns the form-class label associated to the ideal.
    pub fn form_label(&self) -> &IdealFormCorrespondence {
        &self.form_label
    }

    /// Returns how the walk direction was certified.
    pub fn direction_certification(&self) -> CraterDirectionCertification {
        self.direction_certification
    }

    /// Attaches a user-supplied crater orientation witness to this report.
    ///
    /// This produces a separate oriented wrapper instead of changing the
    /// original report in place. The distinction keeps the graph-deterministic
    /// labeled walk available while recording that later interpretation now
    /// depends on an explicit user witness.
    pub fn with_user_orientation(
        self,
        witness: CraterOrientationWitness,
    ) -> Result<OrientedLabeledCraterWalkReport, CraterOrientationWitnessError> {
        OrientedLabeledCraterWalkReport::new(self, witness)
    }
}

/// Failure modes for building a labeled crater walk from an isogeny graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LabeledCraterWalkError {
    /// The graph could not produce the requested local crater report.
    CraterSearch(VolcanoSearchError),
    /// The supplied ideal and class group did not match the crater label data.
    Label(CraterIdealLabelError),
}

impl fmt::Display for LabeledCraterWalkError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CraterSearch(error) => {
                write!(formatter, "could not build crater report: {error}")
            }
            Self::Label(error) => write!(formatter, "could not label crater walk: {error}"),
        }
    }
}

impl std::error::Error for LabeledCraterWalkError {}

impl From<VolcanoSearchError> for LabeledCraterWalkError {
    fn from(error: VolcanoSearchError) -> Self {
        Self::CraterSearch(error)
    }
}

impl From<CraterIdealLabelError> for LabeledCraterWalkError {
    fn from(error: CraterIdealLabelError) -> Self {
        Self::Label(error)
    }
}
