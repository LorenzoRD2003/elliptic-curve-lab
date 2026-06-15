use crate::elliptic_curves::AffinePoint;
use crate::fields::traits::{EnumerableFiniteField, SqrtField};

/// Comparison summary between division-polynomial-based torsion recovery and
/// direct small-group enumeration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TorsionComparisonReport<F: EnumerableFiniteField + SqrtField> {
    n: usize,
    polynomial_candidate_count: usize,
    polynomial_n_torsion_count: usize,
    enumerated_n_torsion_count: usize,
    exact_order_polynomial_count: usize,
    exact_order_enumerated_count: usize,
    missing_from_polynomial: Vec<AffinePoint<F>>,
    extra_from_polynomial: Vec<AffinePoint<F>>,
}

impl<F: EnumerableFiniteField + SqrtField> TorsionComparisonReport<F> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        n: usize,
        polynomial_candidate_count: usize,
        polynomial_n_torsion_count: usize,
        enumerated_n_torsion_count: usize,
        exact_order_polynomial_count: usize,
        exact_order_enumerated_count: usize,
        missing_from_polynomial: Vec<AffinePoint<F>>,
        extra_from_polynomial: Vec<AffinePoint<F>>,
    ) -> Self {
        Self {
            n,
            polynomial_candidate_count,
            polynomial_n_torsion_count,
            enumerated_n_torsion_count,
            exact_order_polynomial_count,
            exact_order_enumerated_count,
            missing_from_polynomial,
            extra_from_polynomial,
        }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    /// Affine points returned by the raw division-polynomial candidate pipeline.
    pub fn polynomial_candidate_count(&self) -> usize {
        self.polynomial_candidate_count
    }

    /// Points kept after the public torsion-validation layer.
    pub fn polynomial_n_torsion_count(&self) -> usize {
        self.polynomial_n_torsion_count
    }

    /// Non-identity affine points satisfying `[n]P = O` by direct group traversal.
    pub fn enumerated_n_torsion_count(&self) -> usize {
        self.enumerated_n_torsion_count
    }

    /// Exact-order points recovered from the division-polynomial pipeline.
    pub fn exact_order_polynomial_count(&self) -> usize {
        self.exact_order_polynomial_count
    }

    /// Exact-order points found by exhaustive group enumeration.
    pub fn exact_order_enumerated_count(&self) -> usize {
        self.exact_order_enumerated_count
    }

    pub fn missing_from_polynomial(&self) -> &[AffinePoint<F>] {
        &self.missing_from_polynomial
    }

    pub fn extra_from_polynomial(&self) -> &[AffinePoint<F>] {
        &self.extra_from_polynomial
    }
}
