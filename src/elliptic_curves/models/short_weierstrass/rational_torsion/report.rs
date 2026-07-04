use num_rational::BigRational;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::rational_torsion::{
        RationalTorsionError, RationalTorsionGroup, enumeration::LutzNagellCandidateReport,
        integral_model::RationalIntegralModel,
        reduction_mod_p::rational_points_from_integral_model, strategy::RationalTorsionStrategy,
        verification::VerifiedRationalTorsion,
    },
};
use crate::fields::Q;

/// Educational report for a completed `E(Q)_tors` computation.
///
/// The point list is the canonical payload. Summary data such as the group
/// classification and candidate counts describe how that list was obtained,
/// but future implementation passes should continue deriving user-facing
/// quantities from `points` whenever possible.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RationalTorsionReport {
    original_curve: ShortWeierstrassCurve<Q>,
    integral_model: ShortWeierstrassCurve<Q>,
    scale: BigRational,
    strategy: RationalTorsionStrategy,
    group: RationalTorsionGroup,
    points: Vec<AffinePoint<Q>>,
    lutz_nagell_candidate_count: Option<usize>,
}

impl RationalTorsionReport {
    /// Builds a completed Lutz-Nagell rational-torsion report from an integral
    /// model.
    ///
    /// This composes the current exact route:
    ///
    /// 1. enumerate Lutz-Nagell candidates on the integral companion;
    /// 2. keep exactly the candidates killed by a possible Mazur point order;
    /// 3. classify the verified finite group;
    /// 4. transport the verified points back to the original source curve.
    pub(crate) fn from_lutz_nagell_integral_model(
        model: &RationalIntegralModel,
    ) -> Result<Self, RationalTorsionError> {
        let candidates = LutzNagellCandidateReport::from_integral_model(model)?;
        let verified = VerifiedRationalTorsion::from_candidates(model.curve(), &candidates)?;
        let source_points = verified
            .points()
            .iter()
            .map(|point| model.to_source_point(point))
            .collect::<Result<Vec<_>, _>>()?;

        Self::new(
            model.source_curve().clone(),
            model.curve().clone(),
            model.scale().clone(),
            RationalTorsionStrategy::LutzNagell,
            verified.group(),
            source_points,
            Some(verified.candidate_count()),
        )
    }

    /// Builds a completed good-reduction/Hensel rational-torsion report from
    /// an integral model.
    pub(crate) fn from_good_reduction_hensel_integral_model(
        model: &RationalIntegralModel,
    ) -> Result<Self, RationalTorsionError> {
        let (group, lifted_points) = rational_points_from_integral_model(model)?;
        let source_points = lifted_points
            .iter()
            .map(|point| model.to_source_point(point))
            .collect::<Result<Vec<_>, _>>()?;

        Self::new(
            model.source_curve().clone(),
            model.curve().clone(),
            model.scale().clone(),
            RationalTorsionStrategy::GoodReductionHensel,
            group,
            source_points,
            None,
        )
    }

    /// Builds a rational-torsion report from already-certified data.
    ///
    /// The constructor checks the basic accounting invariants: the point list
    /// must have the cardinality of the classified group. For the Lutz-Nagell
    /// strategy, the checked candidate count cannot be smaller than the
    /// accepted point count.
    pub(crate) fn new(
        original_curve: ShortWeierstrassCurve<Q>,
        integral_model: ShortWeierstrassCurve<Q>,
        scale: BigRational,
        strategy: RationalTorsionStrategy,
        group: RationalTorsionGroup,
        points: Vec<AffinePoint<Q>>,
        lutz_nagell_candidate_count: Option<usize>,
    ) -> Result<Self, RationalTorsionError> {
        let point_count = points.len();
        let group_cardinality = group.cardinality();
        if point_count != group_cardinality {
            return Err(RationalTorsionError::InconsistentReportGroup {
                group_cardinality,
                point_count,
            });
        }
        if let Some(candidate_count) = lutz_nagell_candidate_count
            && candidate_count < point_count
        {
            return Err(RationalTorsionError::InvalidCandidateAccounting {
                candidate_count,
                point_count,
            });
        }

        Ok(Self {
            original_curve,
            integral_model,
            scale,
            strategy,
            group,
            points,
            lutz_nagell_candidate_count,
        })
    }

    /// Returns the input curve whose torsion subgroup was classified.
    pub fn original_curve(&self) -> &ShortWeierstrassCurve<Q> {
        &self.original_curve
    }

    /// Returns the integral companion model used by the selected strategy.
    pub fn integral_model(&self) -> &ShortWeierstrassCurve<Q> {
        &self.integral_model
    }

    /// Returns the scaling factor `u` for the integral-model transport.
    pub fn scale(&self) -> &BigRational {
        &self.scale
    }

    /// Returns the strategy used to compute this report.
    pub fn strategy(&self) -> RationalTorsionStrategy {
        self.strategy
    }

    /// Returns the Mazur-shape classification of `E(Q)_tors`.
    pub fn group(&self) -> RationalTorsionGroup {
        self.group
    }

    /// Returns the certified rational torsion points, including `O`.
    pub fn points(&self) -> &[AffinePoint<Q>] {
        &self.points
    }

    /// Returns how many Lutz-Nagell candidates were checked, when this report
    /// was produced by the Lutz-Nagell strategy.
    pub fn lutz_nagell_candidate_count(&self) -> Option<usize> {
        self.lutz_nagell_candidate_count
    }

    /// Returns how many Lutz-Nagell candidates were rejected as non-torsion,
    /// when this report was produced by the Lutz-Nagell strategy.
    pub fn lutz_nagell_rejected_candidate_count(&self) -> Option<usize> {
        self.lutz_nagell_candidate_count
            .map(|candidate_count| candidate_count - self.points.len())
    }
}

impl ShortWeierstrassCurve<Q> {
    /// Computes and classifies the rational torsion subgroup `E(ℚ)_tors`.
    ///
    /// The strategy is explicit so callers can choose between the classical
    /// Lutz-Nagell route and the good-reduction/Hensel route.
    pub fn rational_torsion_by(
        &self,
        strategy: RationalTorsionStrategy,
    ) -> Result<RationalTorsionReport, RationalTorsionError> {
        let model = RationalIntegralModel::from_curve(self.clone())?;
        match strategy {
            RationalTorsionStrategy::LutzNagell => {
                RationalTorsionReport::from_lutz_nagell_integral_model(&model)
            }
            RationalTorsionStrategy::GoodReductionHensel => {
                RationalTorsionReport::from_good_reduction_hensel_integral_model(&model)
            }
        }
    }
}
