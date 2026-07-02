use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve,
    traits::{GroupCurveModel, HasProjectiveModel, ProjectiveGroupCurveModel, ScalarInput},
};
use crate::fields::traits::*;
use num_bigint::{BigInt, Sign};

impl<F: Field> GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    // The affine public API now delegates to the native homogeneous projective
    // core so the family has one canonical arithmetic engine. The actual
    // formulas used by that engine are documented on the projective helpers.

    fn negated_y_coordinate(&self, x: &F::Elem, y: &F::Elem) -> F::Elem {
        F::sub(&F::neg(y), &F::add(&F::mul(self.a1(), x), self.a3()))
    }

    /// Adds two already-validated affine points by lifting them to the native
    /// projective surface, applying the homogeneous group law there, and
    /// recovering the affine result.
    ///
    /// Complexity: `O(1)` projective group-law work plus two affine lifts and
    /// one affine recovery.
    fn add_unchecked(
        &self,
        left: &AffinePoint<F>,
        right: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        let left_projective = self.to_projective(left)?;
        let right_projective = self.to_projective(right)?;
        let sum_projective = self.add_projective(&left_projective, &right_projective)?;
        self.to_affine_projective(&sum_projective)
    }

    /// Doubles one already-validated affine point through the native
    /// homogeneous projective core.
    ///
    /// Complexity: `O(1)` projective doubling work plus one affine lift and one
    /// affine recovery.
    fn double_unchecked(&self, point: &AffinePoint<F>) -> Result<AffinePoint<F>, CurveError> {
        let projective = self.to_projective(point)?;
        let doubled_projective = self.double_projective(&projective)?;
        self.to_affine_projective(&doubled_projective)
    }

    /// Multiplies one already-validated affine point by a non-negative scalar
    /// through the native projective add-and-double core.
    ///
    /// Complexity: `Θ(log n)` projective additions/doublings for scalar `n`,
    /// plus one affine lift and one affine recovery.
    fn mul_scalar_unchecked(
        &self,
        point: &AffinePoint<F>,
        scalar: impl ScalarInput,
    ) -> Result<AffinePoint<F>, CurveError> {
        let projective = self.to_projective(point)?;
        let multiple_projective = self.mul_scalar_projective(&projective, scalar)?;
        self.to_affine_projective(&multiple_projective)
    }
}

impl<F: Field> GroupCurveModel for GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn neg(&self, point: &Self::Point) -> Self::Point {
        match point {
            AffinePoint::Infinity => AffinePoint::Infinity,
            AffinePoint::Finite { x, y } => {
                let negated = AffinePoint::new(x.clone(), self.negated_y_coordinate(x, y));
                debug_assert!(
                    !self.contains_affine_point(point) || self.contains_affine_point(&negated),
                    "general-Weierstrass negation should preserve on-curve inputs"
                );
                negated
            }
        }
    }

    fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(left) || !self.contains_affine_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.add_unchecked(left, right)
    }

    fn sub(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(left) || !self.contains_affine_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let negated = self.neg(right);
        self.add(left, &negated)
    }

    fn double(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.double_unchecked(point)
    }

    fn mul_scalar(
        &self,
        point: &Self::Point,
        scalar: impl ScalarInput,
    ) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.mul_scalar_unchecked(point, scalar)
    }

    fn mul_scalar_signed(
        &self,
        point: &Self::Point,
        scalar: impl Into<BigInt>,
    ) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let scalar = scalar.into();
        if scalar.sign() == Sign::Minus {
            let negated = self.neg(point);
            self.mul_scalar_unchecked(
                &negated,
                (-scalar)
                    .to_biguint()
                    .expect("negated negative scalar should be non-negative"),
            )
        } else {
            self.mul_scalar_unchecked(
                point,
                scalar
                    .to_biguint()
                    .expect("non-negative scalar should convert to BigUint"),
            )
        }
    }
}
