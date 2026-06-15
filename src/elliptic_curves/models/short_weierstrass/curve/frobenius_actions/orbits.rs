use std::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    frobenius::orbit::{
        FrobeniusOrbit, orbit_from_successor_by_key, partition_point_orbits_by_key,
    },
    traits::{CurveModel, EnumerableCurveModel},
};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Applies the absolute Frobenius `π_p^k` to a point on this curve.
    ///
    /// If `P = (x, y)` lies on `E`, the returned point is
    /// `(x^(p^k), y^(p^k))`, which lies on the twist
    /// `self.frobenius_twist_power(power)`.
    pub fn absolute_frobenius_power_point(
        &self,
        point: &AffinePoint<F>,
        power: u32,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        Ok(point
            .map_coordinates::<F, _>(|coordinate| Self::frobenius_power_element(coordinate, power)))
    }

    /// Computes the orbit of one rational point under the absolute
    /// Frobenius `π_p^k`.
    pub fn absolute_frobenius_orbit(
        &self,
        point: &AffinePoint<F>,
        power: u32,
    ) -> Result<FrobeniusOrbit<AffinePoint<F>>, CurveError>
    where
        F::Elem: Hash,
    {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        self.absolute_frobenius_preserves_curve(power)?;

        let bound = self.absolute_frobenius_period_bound(power);
        orbit_from_successor_by_key(
            point.clone(),
            bound as usize,
            |orbit_point| orbit_point.clone(),
            |current| self.absolute_frobenius_power_point(current, power),
        )
    }

    /// Partitions the represented rational point set into orbits of the
    /// absolute Frobenius `π_p^k`.
    pub fn absolute_frobenius_orbits_on_points(
        &self,
        power: u32,
    ) -> Result<Vec<FrobeniusOrbit<AffinePoint<F>>>, CurveError>
    where
        F: EnumerableFiniteField + SqrtField,
        F::Elem: Hash,
    {
        self.absolute_frobenius_preserves_curve(power)?;

        partition_point_orbits_by_key(
            self.points(),
            |orbit_point| orbit_point.clone(),
            |point| self.absolute_frobenius_orbit(point, power),
        )
    }
}
