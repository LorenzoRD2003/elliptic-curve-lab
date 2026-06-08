use crate::elliptic_curves::frobenius::FrobeniusTrace;
use crate::elliptic_curves::{
    CurveError, FrobeniusTraceCurveModel, ShortWeierstrassQuadraticTwist,
};
use crate::fields::{EnumerableFiniteField, FiniteField, SqrtField};

/// Frobenius relation between a curve and a chosen quadratic twist over `F_q`.
///
/// For a non-trivial quadratic twist `E'` of `E` over `F_q`, the classical
/// point-count identity is `#E(F_q) + #E'(F_q) = 2q + 2`.
///
/// Equivalently, if `#E(F_q) = q + 1 - t`, then `#E'(F_q) = q + 1 + t`,
/// so the two Frobenius traces satisfy `t' = -t`.
///
/// This report stores the two Frobenius-trace packages as the primary data and
/// records whether the expected order relation holds for the chosen twist.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuadraticTwistFrobeniusRelation {
    original: FrobeniusTrace,
    twist: FrobeniusTrace,
    sum_orders: u128,
    expected_sum: u128,
    holds: bool,
}

impl QuadraticTwistFrobeniusRelation {
    /// Returns the Frobenius trace package of the original curve `E`.
    pub fn original(&self) -> &FrobeniusTrace {
        &self.original
    }

    /// Returns the Frobenius trace package of the chosen twist `E'`.
    pub fn twist(&self) -> &FrobeniusTrace {
        &self.twist
    }

    /// Returns `#E(F_q) + #E'(F_q)`.
    pub fn sum_orders(&self) -> u128 {
        self.sum_orders
    }

    /// Returns the expected value `2q + 2`.
    pub fn expected_sum(&self) -> u128 {
        self.expected_sum
    }

    /// Returns whether `#E(F_q) + #E'(F_q) = 2q + 2`.
    pub fn holds(&self) -> bool {
        self.holds
    }

    /// Returns whether the two Frobenius traces satisfy `t' = -t`.
    pub fn trace_negation_holds(&self) -> bool {
        self.twist.trace() == -self.original.trace()
    }
}

impl<F: EnumerableFiniteField + SqrtField + FiniteField> ShortWeierstrassQuadraticTwist<F> {
    /// Computes the Frobenius relation between `E` and the stored twist `E'`.
    ///
    /// If the chosen twist factor is genuinely quadratic over `F_q`, one
    /// expects `#E(F_q) + #E'(F_q) = 2q + 2`, and equivalently `t' = -t`.
    ///
    /// The current implementation derives both traces from exhaustive point
    /// counts on the two curves and then compares the resulting invariants.
    ///
    /// Complexity: `Θ(1)`
    pub fn frobenius_relation(&self) -> Result<QuadraticTwistFrobeniusRelation, CurveError> {
        let original = self.original().frobenius_trace()?;
        let twist = self.twist().frobenius_trace()?;

        let sum_orders = u128::from(original.curve_order()) + u128::from(twist.curve_order());
        let field_order = original.field_order();
        let expected_sum = field_order
            .checked_mul(2)
            .and_then(|double| double.checked_add(2))
            .expect("enumerable finite-field Frobenius relation should keep 2q + 2 inside u128");
        let holds = sum_orders == expected_sum;

        Ok(QuadraticTwistFrobeniusRelation {
            original,
            twist,
            sum_orders,
            expected_sum,
            holds,
        })
    }
}
