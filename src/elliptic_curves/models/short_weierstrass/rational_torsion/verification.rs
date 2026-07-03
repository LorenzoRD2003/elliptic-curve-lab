use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::rational_torsion::{
        RationalTorsionError, RationalTorsionGroup,
        enumeration::{LutzNagellCandidateReport, sort_affine_points},
        mazur::MAZUR_CYCLIC_ORDERS,
    },
    traits::{CurveModel, GroupCurveModel},
};
use crate::fields::Q;

#[derive(Clone, Debug, PartialEq, Eq)]
struct VerifiedTorsionPoint {
    point: AffinePoint<Q>,
    order: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct VerifiedRationalTorsion {
    group: RationalTorsionGroup,
    points: Vec<AffinePoint<Q>>,
    candidate_count: usize,
}

impl VerifiedRationalTorsion {
    /// Verifies and classifies the rational torsion points among Lutz-Nagell
    /// candidates.
    ///
    /// A candidate survives exactly when it is `O`, or when it is killed by
    /// one of Mazur's possible non-identity rational point orders:
    /// `2, ..., 10, 12`. This avoids multiplying non-torsion candidates by the
    /// much larger exponent `27720` over `Q`, where coefficient growth can be
    /// severe.
    ///
    /// Complexity: `Θ(C + T log T)`, where `C` is the number of candidates and
    /// `T` is the number of surviving torsion points.
    pub(crate) fn from_candidates(
        curve: &ShortWeierstrassCurve<Q>,
        candidates: &LutzNagellCandidateReport,
    ) -> Result<Self, RationalTorsionError> {
        let mut verified_points = Vec::new();
        for candidate in candidates.points() {
            if let Some(order) = curve.exact_mazur_order(candidate)? {
                verified_points.push(VerifiedTorsionPoint {
                    point: candidate.clone(),
                    order,
                });
            }
        }

        let (mut torsion_points, point_orders): (Vec<_>, Vec<_>) = verified_points
            .into_iter()
            .map(|verified| (verified.point, verified.order))
            .unzip();
        sort_affine_points(&mut torsion_points);
        let group = RationalTorsionGroup::from_verified_point_orders(&point_orders)?;

        Ok(Self {
            group,
            points: torsion_points,
            candidate_count: candidates.candidate_count(),
        })
    }

    pub(crate) fn group(&self) -> RationalTorsionGroup {
        self.group
    }

    pub(crate) fn points(&self) -> &[AffinePoint<Q>] {
        &self.points
    }

    pub(crate) fn candidate_count(&self) -> usize {
        self.candidate_count
    }
}

impl ShortWeierstrassCurve<Q> {
    pub(crate) fn exact_mazur_order(
        &self,
        point: &AffinePoint<Q>,
    ) -> Result<Option<usize>, RationalTorsionError> {
        if !self.contains(point) {
            return Err(RationalTorsionError::from(
                crate::elliptic_curves::CurveError::PointNotOnCurve,
            ));
        }
        if self.is_identity(point) {
            return Ok(Some(1));
        }

        for order in MAZUR_CYCLIC_ORDERS {
            if self.is_identity(&self.mul_scalar(point, *order)?) {
                return Ok(Some(*order));
            }
        }
        Ok(None)
    }
}
