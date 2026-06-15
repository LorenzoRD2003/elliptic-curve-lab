use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::{DivisionPolynomialError, TorsionComparisonReport},
    traits::{CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel},
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};

fn difference_of_point_sets<F: EnumerableFiniteField + SqrtField>(
    left: &[AffinePoint<F>],
    right: &[AffinePoint<F>],
) -> Vec<AffinePoint<F>> {
    let mut difference = Vec::new();
    for point in left {
        if !right.iter().any(|candidate| candidate == point) {
            difference.push(point.clone());
        }
    }
    difference
}

impl<F: EnumerableFiniteField + SqrtField> ShortWeierstrassCurve<F> {
    /// Compares division-polynomial torsion recovery against direct exhaustive
    /// group enumeration.
    ///
    /// This helper is meant for the current small finite educational setting where
    /// the full rational point set can be enumerated honestly. It summarizes:
    ///
    /// - how many affine points arise as raw division-polynomial candidates
    /// - how many survive the public torsion validation layer
    /// - how many are found by direct `[n]P = O` enumeration
    /// - how the exact-order-`n` sets compare
    ///
    /// The `missing_from_polynomial` and `extra_from_polynomial` fields compare
    /// the exact-order sets, since that is the most informative notion of
    /// correctness once candidates and lower-order contamination are taken into
    /// account.
    pub fn compare_division_polynomial_torsion_with_enumeration(
        &self,
        n: usize,
    ) -> Result<TorsionComparisonReport<F>, DivisionPolynomialError> {
        let polynomial_candidates = self.torsion_candidates_from_division_polynomial(n)?;
        let polynomial_n_torsion = self.torsion_points_from_division_polynomial(n)?;
        let exact_order_polynomial = self.exact_n_torsion_points_from_division_polynomial(n)?;

        let enumerated_n_torsion: Vec<_> = self
            .points()
            .into_iter()
            .filter(|point| !self.is_identity(point) && self.is_torsion_point(point, n as u64))
            .collect();

        let exact_order_enumerated = self
            .points_of_exact_order(n)
            .map_err(DivisionPolynomialError::Curve)?;

        let missing_from_polynomial =
            difference_of_point_sets(&exact_order_enumerated, &exact_order_polynomial);
        let extra_from_polynomial =
            difference_of_point_sets(&exact_order_polynomial, &exact_order_enumerated);

        Ok(TorsionComparisonReport::new(
            n,
            polynomial_candidates.len(),
            polynomial_n_torsion.len(),
            enumerated_n_torsion.len(),
            exact_order_polynomial.len(),
            exact_order_enumerated.len(),
            missing_from_polynomial,
            extra_from_polynomial,
        ))
    }
}
