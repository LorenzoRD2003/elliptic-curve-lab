use crate::elliptic_curves::{
    AffinePoint, CurveError, TwistedEdwardsCurve,
    traits::{GroupCurveModel, HasProjectiveModel, ProjectiveGroupCurveModel},
};
use crate::fields::traits::Field;

impl<F: Field> TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    fn affine_neg_x(&self, x: &F::Elem) -> F::Elem {
        F::neg(x)
    }
}

impl<F: Field> GroupCurveModel for TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    fn neg(&self, point: &Self::Point) -> Self::Point {
        match point {
            AffinePoint::Infinity => AffinePoint::Infinity,
            AffinePoint::Finite { x, y } => {
                let negated = AffinePoint::new(self.affine_neg_x(x), y.clone());
                debug_assert!(
                    !self.contains_affine_point(point) || self.contains_affine_point(&negated),
                    "twisted-Edwards negation should preserve on-curve inputs"
                );
                negated
            }
        }
    }

    fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(left) || !self.contains_affine_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        // The executable core now lives in the extended-coordinate projective
        // layer. If the native projective sum leaves the affine chart and
        // returns one valid point with Z = 0, the group law itself has still
        // succeeded; the only failing step is the final recovery back to
        // affine coordinates.
        let left_projective = self.to_projective(left)?;
        let right_projective = self.to_projective(right)?;
        let sum = self.add_projective(&left_projective, &right_projective)?;
        self.to_affine_projective(&sum)
    }

    fn sub(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(left) || !self.contains_affine_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let left_projective = self.to_projective(left)?;
        let right_projective = self.to_projective(right)?;
        let difference =
            self.add_projective(&left_projective, &self.neg_projective(&right_projective))?;
        self.to_affine_projective(&difference)
    }

    fn double(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        // As in add(), a valid projective double may leave the affine chart.
        // In that case the honest error comes from affine recovery, not from
        // the native projective doubling formula itself.
        let projective = self.to_projective(point)?;
        let doubled = self.double_projective(&projective)?;
        self.to_affine_projective(&doubled)
    }

    fn mul_scalar(&self, point: &Self::Point, scalar: u64) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let projective = self.to_projective(point)?;
        let multiple = self.mul_scalar_projective(&projective, scalar)?;
        self.to_affine_projective(&multiple)
    }

    fn mul_scalar_signed(
        &self,
        point: &Self::Point,
        scalar: i64,
    ) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let projective = self.to_projective(point)?;
        if scalar < 0 {
            let multiple = self
                .mul_scalar_projective(&self.neg_projective(&projective), scalar.unsigned_abs())?;
            self.to_affine_projective(&multiple)
        } else {
            let multiple = self.mul_scalar_projective(&projective, scalar as u64)?;
            self.to_affine_projective(&multiple)
        }
    }
}
