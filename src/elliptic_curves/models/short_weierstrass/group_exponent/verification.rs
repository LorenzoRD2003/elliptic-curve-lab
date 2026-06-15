use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    frobenius::group_order::{GroupOrderReport, GroupOrderStrategy},
    short_weierstrass::group_exponent::ExponentAccumulationReport,
    traits::CurveModel,
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};

/// Group-order-side verification of one accumulated exponent lower bound.
///
/// This report does not certify the exponent itself. It records whether the
/// Hasse interval attached to one chosen group-order route contains a unique
/// multiple of the supplied lower bound, which would force one group order
/// `#E(F_q)` compatible with that lower bound.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExponentLowerBoundGroupOrderVerification {
    exponent_lower_bound: BigUint,
    group_order_report: GroupOrderReport,
}

impl ExponentLowerBoundGroupOrderVerification {
    pub(crate) fn new(exponent_lower_bound: BigUint, group_order_report: GroupOrderReport) -> Self {
        Self {
            exponent_lower_bound,
            group_order_report,
        }
    }

    /// Returns the lower bound being checked.
    pub fn exponent_lower_bound(&self) -> &BigUint {
        &self.exponent_lower_bound
    }

    /// Returns the group-order report that supplied the Hasse interval.
    pub fn group_order_report(&self) -> &GroupOrderReport {
        &self.group_order_report
    }

    fn unique_group_order_multiple_in_hasse_interval(&self) -> Option<u128> {
        self.exponent_lower_bound.to_u128().and_then(|lower_bound| {
            self.group_order_report
                .hasse_interval()
                .unique_multiple_of(lower_bound)
        })
    }

    /// Returns the verified group order when the Hasse interval forces one
    /// unique multiple of the lower bound.
    pub fn verified_group_order(&self) -> Option<u128> {
        self.unique_group_order_multiple_in_hasse_interval()
    }
}

impl<F: FiniteField + EnumerableFiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    pub(crate) fn verify_exponent_lower_bound_by_group_order_report(
        &self,
        exponent_lower_bound: BigUint,
        group_order_report: GroupOrderReport,
    ) -> ExponentLowerBoundGroupOrderVerification {
        ExponentLowerBoundGroupOrderVerification::new(exponent_lower_bound, group_order_report)
    }

    /// Verifies one accumulated exponent lower bound against a chosen
    /// group-order route.
    ///
    /// This method is intentionally separate from [`Self::group_exponent_by`]:
    /// the random-point exponent route stays a pure lower-bound accumulator,
    /// while this helper uses one explicit [`GroupOrderStrategy`] to ask
    /// whether the resulting Hasse interval `H(q)` contains a unique multiple
    /// of that lower bound.
    ///
    /// If the returned report has `verified_group_order = Some(N)`, then the
    /// Hasse interval for the chosen group-order route contains exactly one
    /// multiple of the lower bound, namely `N`. This certifies one possible
    /// group order `#E(F_q)`, not the exponent itself.
    ///
    /// The intended workflow is to pass an
    /// [`ExponentAccumulationReport<AffinePoint<F>>`] produced from this same
    /// curve. The method rejects obviously incompatible reports whose sampled
    /// points do not lie on the current curve.
    pub fn verify_exponent_lower_bound_by_group_order(
        &self,
        accumulation: &ExponentAccumulationReport<AffinePoint<F>>,
        strategy: GroupOrderStrategy,
    ) -> Result<ExponentLowerBoundGroupOrderVerification, CurveError> {
        for step in accumulation.steps() {
            if !self.contains(step.point()) {
                return Err(CurveError::PointNotOnCurve);
            }
        }

        Ok(self.verify_exponent_lower_bound_by_group_order_report(
            accumulation.exponent_lower_bound().clone(),
            self.group_order_by(strategy)?,
        ))
    }
}
