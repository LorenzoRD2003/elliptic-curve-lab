use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::DivisionPolynomialError,
    traits::{FiniteGroupCurveModel, GroupCurveModel},
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};

impl<F: EnumerableFiniteField + SqrtField> ShortWeierstrassCurve<F> {
    /// Returns the rational affine points kept after a light torsion validation
    /// pass on top of the division-polynomial candidates.
    ///
    /// For odd `n`, the current search already lands exactly on the affine
    /// `n`-torsion candidates detected by `ψ_n(x)`, so this step is just the
    /// identity on the candidate set.
    ///
    /// For even `n`, the `x`-criterion only sees the stripped factor in
    /// `ψ_n = y ε_n(x)`, so affine points with `y = 0` may appear for reasons
    /// unrelated to exact `n`-torsion. The current educational validation keeps
    /// nonzero-`y` points directly and checks `[n]P = O` for the zero-`y` cases.
    ///
    /// Complexity: `Θ(q²)` affine candidate work plus, in the even case,
    /// additional torsion checks on the `y = 0` candidates using the current
    /// small-group group-law backend.
    pub fn torsion_points_from_division_polynomial(
        &self,
        n: usize,
    ) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
        let candidates = self.torsion_candidates_from_division_polynomial(n)?;
        if !self.division_polynomial_uses_y_factor(n)? {
            return Ok(candidates);
        }

        let mut points = Vec::new();
        for point in candidates {
            if !self.point_has_zero_y(&point) || self.is_torsion_point(&point, n as u64) {
                points.push(point);
            }
        }
        Ok(points)
    }

    /// Returns the rational affine points of exact order `n` recovered from the
    /// division-polynomial search pipeline.
    ///
    /// This route first builds the public torsion-point candidate set and then
    /// filters by exact order through the ambient small-group model.
    ///
    /// Complexity: `Θ(q²)` candidate work plus one exact-order check per kept
    /// affine candidate.
    pub fn exact_n_torsion_points_from_division_polynomial(
        &self,
        n: usize,
    ) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
        let candidates = self.torsion_points_from_division_polynomial(n)?;
        let mut points = Vec::new();

        for point in candidates {
            if self
                .point_has_exact_order(&point, n)
                .map_err(|error| match error {
                    CurveError::InvalidTorsionOrder { order: _ } => {
                        DivisionPolynomialError::ZeroIndex
                    }
                    other => DivisionPolynomialError::Curve(other),
                })?
            {
                points.push(point);
            }
        }

        Ok(points)
    }
}
