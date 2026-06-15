use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, EllipticFunctionTruncation,
    zeta::value::WeierstrassZetaApprox,
};
use crate::fields::complex_approx::ComplexApprox;

impl ComplexLattice {
    /// Approximates the Weierstrass `ζ`-function attached to `self`.
    pub fn weierstrass_zeta(
        &self,
        z: Complex64,
        truncation: EllipticFunctionTruncation,
    ) -> Result<WeierstrassZetaApprox, AnalyticCurveError> {
        let canonical_coordinate = self.reduce_complex_point_to_fundamental_coordinates(z)?;
        let canonical_z = self.point_from_fundamental_coordinates(canonical_coordinate);
        let tolerance = ComplexApprox::default_tolerance();

        if ComplexApprox::is_zero_with_tolerance(&canonical_z, tolerance) {
            return Err(AnalyticCurveError::PointTooCloseToLatticePoint);
        }

        let lattice_points = self.nonzero_lattice_points_in_box(truncation.radius());
        let one = Complex64::new(1.0, 0.0);
        let mut value = one / z;
        let mut pole_distance = canonical_z.norm();

        for point in &lattice_points {
            let shifted_for_value = z - point.value;
            let shifted_for_pole = canonical_z - point.value;

            if ComplexApprox::is_zero_with_tolerance(&shifted_for_pole, tolerance) {
                return Err(AnalyticCurveError::PointTooCloseToLatticePoint);
            }

            pole_distance = pole_distance.min(shifted_for_pole.norm());
            value += (one / shifted_for_value) + (one / point.value) + (z / point.value.powu(2));
        }

        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(AnalyticCurveError::NumericalComparisonFailed);
        }

        Ok(WeierstrassZetaApprox::new(
            z,
            value,
            truncation,
            lattice_points.len(),
            pole_distance,
        ))
    }
}
