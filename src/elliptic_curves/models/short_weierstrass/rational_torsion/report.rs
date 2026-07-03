use num_rational::BigRational;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve, short_weierstrass::rational_torsion::RationalTorsionGroup,
};
use crate::fields::Q;

/// Educational report for a completed `E(Q)_tors` computation.
///
/// The point list is the canonical payload. Summary data such as the group
/// classification and candidate counts describe how that list was obtained,
/// but future implementation passes should continue deriving user-facing
/// quantities from `points` whenever possible.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RationalTorsionReport {
    original_curve: ShortWeierstrassCurve<Q>,
    integral_model: ShortWeierstrassCurve<Q>,
    scale: BigRational,
    group: RationalTorsionGroup,
    points: Vec<AffinePoint<Q>>,
    candidate_count: usize,
    rejected_candidate_count: usize,
}

impl RationalTorsionReport {
    /// Builds a rational-torsion report from already-certified data.
    pub(crate) fn new(
        original_curve: ShortWeierstrassCurve<Q>,
        integral_model: ShortWeierstrassCurve<Q>,
        scale: BigRational,
        group: RationalTorsionGroup,
        points: Vec<AffinePoint<Q>>,
        candidate_count: usize,
        rejected_candidate_count: usize,
    ) -> Self {
        Self {
            original_curve,
            integral_model,
            scale,
            group,
            points,
            candidate_count,
            rejected_candidate_count,
        }
    }

    /// Returns the input curve whose torsion subgroup was classified.
    pub(crate) fn original_curve(&self) -> &ShortWeierstrassCurve<Q> {
        &self.original_curve
    }

    /// Returns the integral companion model used for Lutz-Nagell search.
    pub(crate) fn integral_model(&self) -> &ShortWeierstrassCurve<Q> {
        &self.integral_model
    }

    /// Returns the scaling factor `u` for the integral-model transport.
    pub(crate) fn scale(&self) -> &BigRational {
        &self.scale
    }

    /// Returns the Mazur-shape classification of `E(Q)_tors`.
    pub(crate) fn group(&self) -> RationalTorsionGroup {
        self.group
    }

    /// Returns the certified rational torsion points, including `O`.
    pub(crate) fn points(&self) -> &[AffinePoint<Q>] {
        &self.points
    }

    /// Returns how many integral candidates were checked.
    pub(crate) fn candidate_count(&self) -> usize {
        self.candidate_count
    }

    /// Returns how many checked candidates were rejected as non-torsion.
    pub(crate) fn rejected_candidate_count(&self) -> usize {
        self.rejected_candidate_count
    }
}
