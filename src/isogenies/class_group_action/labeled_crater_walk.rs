use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::QuadraticClassGroup,
    quadratic_ideals::{IdealFormCorrespondence, PrimeNormIdeal},
};
use crate::isogenies::{
    class_group_action::{CraterIdealLabelError, CraterIdealLabelReport, CraterWalkReport},
    graphs::{IsogenyGraphNodeId, endomorphisms::CraterReport},
};

/// Certification status for the direction used by a labeled crater walk.
///
/// The first labeled-walk milestone deliberately records only that the walk
/// direction was chosen deterministically from graph data. It does not certify
/// which horizontal direction corresponds to `𝔭` rather than `\bar{𝔭}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CraterDirectionCertification {
    /// The walk follows the deterministic graph-local rule from [`CraterWalkReport`].
    GraphDeterministic,
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
pub(crate) struct LabeledCraterWalkReport {
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
    pub(crate) fn local_label(&self) -> &CraterIdealLabelReport {
        &self.local_label
    }

    /// Returns the deterministic crater walk.
    pub(crate) fn walk(&self) -> &CraterWalkReport {
        &self.walk
    }

    /// Returns the form-class label associated to the ideal.
    pub(crate) fn form_label(&self) -> &IdealFormCorrespondence {
        &self.form_label
    }

    /// Returns how the walk direction was certified.
    pub(crate) fn direction_certification(&self) -> CraterDirectionCertification {
        self.direction_certification
    }
}
