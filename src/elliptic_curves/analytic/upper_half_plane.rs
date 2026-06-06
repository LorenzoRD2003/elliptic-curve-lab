use num_complex::Complex64;

use crate::fields::ComplexApprox;

use crate::elliptic_curves::analytic::AnalyticCurveError;

/// A point `τ` in the complex upper half-plane
/// `ℍ = {τ ∈ ℂ : Im(τ) > 0}`.
///
/// This type is the analytic-side wrapper around a validated complex number.
/// The stored value is a concrete [`Complex64`], while validation uses the
/// current default tolerance policy from [`ComplexApprox`] to avoid treating
/// numerically noisy boundary values as honest upper-half-plane points.
#[derive(Clone, Debug, PartialEq)]
pub struct UpperHalfPlanePoint {
    tau: Complex64,
}

impl UpperHalfPlanePoint {
    /// Returns the classical point `τ = i`.
    ///
    /// This is the standard fixed point of the inversion
    /// `τ ↦ -1/τ` in the upper half-plane.
    pub fn tau_i() -> UpperHalfPlanePoint {
        Self::new(Complex64::new(0.0, 1.0)).expect("i lies in the upper half-plane")
    }

    /// Returns the classical point `τ = ρ = -1/2 + (√3/2)i`.
    ///
    /// This is the usual hexagonal-lattice parameter and the standard
    /// order-three elliptic point for the modular group.
    pub fn tau_rho() -> UpperHalfPlanePoint {
        Self::new(Complex64::new(-0.5, f64::sqrt(3.0) * 0.5))
            .expect("rho lies in the upper half-plane")
    }

    /// Returns a simple interior sample point for educational experiments.
    ///
    /// Unlike [`Self::tau_i`] and [`Self::tau_rho`], this point is not a
    /// special modular fixed point. It is just a convenient generic example
    /// that stays comfortably away from the real axis.
    pub fn tau_generic_example() -> UpperHalfPlanePoint {
        Self::new(Complex64::new(1.0 / 3.0, 5.0 / 4.0))
            .expect("the generic sample point lies in the upper half-plane")
    }

    /// Builds a validated upper-half-plane point from a complex value.
    ///
    /// Values whose imaginary part is not strictly above the default absolute
    /// tolerance are rejected with [`AnalyticCurveError::TauNotInUpperHalfPlane`].
    pub fn new(tau: Complex64) -> Result<Self, AnalyticCurveError> {
        if Self::is_in_upper_half_plane(&tau) {
            Ok(Self { tau })
        } else {
            Err(AnalyticCurveError::TauNotInUpperHalfPlane)
        }
    }

    /// Builds a validated upper-half-plane point from real and imaginary
    /// coordinates.
    pub fn from_re_im(re: f64, im: f64) -> Result<Self, AnalyticCurveError> {
        Self::new(Complex64::new(re, im))
    }

    /// Returns the stored upper-half-plane parameter `τ`.
    pub fn tau(&self) -> &Complex64 {
        &self.tau
    }

    /// Returns the real part `Re(τ)`.
    pub fn real_part(&self) -> f64 {
        self.tau.re
    }

    /// Returns the imaginary part `Im(τ)`.
    pub fn imaginary_part(&self) -> f64 {
        self.tau.im
    }

    /// Returns the squared complex norm `|τ|^2`.
    pub fn norm_sqr(&self) -> f64 {
        self.tau.norm_sqr()
    }

    /// Returns whether a complex value lies in the numerical upper half-plane.
    ///
    /// The current implementation treats `τ` as inside only when `Im(τ)` is
    /// strictly larger than the default absolute tolerance used by [`ComplexApprox`].
    pub fn is_in_upper_half_plane(tau: &Complex64) -> bool {
        tau.im > ComplexApprox::default_tolerance().absolute
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;
    use proptest::prelude::*;

    use crate::elliptic_curves::analytic::UpperHalfPlanePoint;
    use crate::{elliptic_curves::analytic::AnalyticCurveError, fields::ComplexApprox};

    #[test]
    fn constructor_accepts_values_with_clearly_positive_imaginary_part() {
        let tau = Complex64::new(0.25, 1.5);
        let point = UpperHalfPlanePoint::new(tau).expect("point should be valid");

        assert_eq!(point.tau(), &tau);
    }

    #[test]
    fn constructor_rejects_real_axis_and_lower_half_plane_values() {
        assert_eq!(
            UpperHalfPlanePoint::new(Complex64::new(1.0, 0.0)),
            Err(AnalyticCurveError::TauNotInUpperHalfPlane)
        );
        assert_eq!(
            UpperHalfPlanePoint::new(Complex64::new(-2.0, -0.5)),
            Err(AnalyticCurveError::TauNotInUpperHalfPlane)
        );
    }

    #[test]
    fn constructor_rejects_near_boundary_values_inside_default_tolerance() {
        let tiny_imaginary_part = ComplexApprox::default_tolerance().absolute * 0.25;

        assert_eq!(
            UpperHalfPlanePoint::new(Complex64::new(0.0, tiny_imaginary_part)),
            Err(AnalyticCurveError::TauNotInUpperHalfPlane)
        );
    }

    #[test]
    fn from_re_im_and_coordinate_accessors_are_consistent() {
        let point = UpperHalfPlanePoint::from_re_im(-0.75, 2.0).expect("point should be valid");

        assert_eq!(point.real_part(), -0.75);
        assert_eq!(point.imaginary_part(), 2.0);
        assert_eq!(point.tau(), &Complex64::new(-0.75, 2.0));
    }

    #[test]
    fn norm_sqr_matches_complex_backend_norm() {
        let point = UpperHalfPlanePoint::from_re_im(3.0, 4.0).expect("point should be valid");

        assert_eq!(point.norm_sqr(), 25.0);
    }

    #[test]
    fn membership_helper_uses_the_same_boundary_rule_as_constructor() {
        let tolerance = ComplexApprox::default_tolerance().absolute;

        assert!(!UpperHalfPlanePoint::is_in_upper_half_plane(
            &Complex64::new(0.0, tolerance * 0.5)
        ));
        assert!(UpperHalfPlanePoint::is_in_upper_half_plane(
            &Complex64::new(0.0, tolerance * 2.0)
        ));
    }

    #[test]
    fn tau_i_is_the_standard_imaginary_unit_point() {
        let point = UpperHalfPlanePoint::tau_i();

        assert_eq!(point.tau(), &Complex64::new(0.0, 1.0));
    }

    #[test]
    fn tau_rho_matches_the_hexagonal_lattice_parameter() {
        let point = UpperHalfPlanePoint::tau_rho();

        assert_eq!(point.real_part(), -0.5);
        assert_eq!(point.imaginary_part(), f64::sqrt(3.0) * 0.5);
    }

    #[test]
    fn generic_example_is_valid_and_non_special() {
        let point = UpperHalfPlanePoint::tau_generic_example();

        assert_eq!(point.tau(), &Complex64::new(1.0 / 3.0, 5.0 / 4.0));
        assert_ne!(point, UpperHalfPlanePoint::tau_i());
        assert_ne!(point, UpperHalfPlanePoint::tau_rho());
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]

        #[test]
        fn generic_positive_imaginary_parts_are_accepted(
            re in -5.0f64..5.0,
            im in 0.1f64..5.0,
        ) {
            let point = UpperHalfPlanePoint::from_re_im(re, im).unwrap();

            prop_assert_eq!(point.real_part(), re);
            prop_assert_eq!(point.imaginary_part(), im);
            prop_assert_eq!(point.tau(), &Complex64::new(re, im));
            prop_assert_eq!(point.norm_sqr(), point.tau().norm_sqr());
            prop_assert!(UpperHalfPlanePoint::is_in_upper_half_plane(point.tau()));
        }

        #[test]
        fn generic_non_positive_imaginary_parts_are_rejected(
            re in -5.0f64..5.0,
            im in -5.0f64..=ComplexApprox::default_tolerance().absolute,
        ) {
            prop_assert_eq!(
                UpperHalfPlanePoint::from_re_im(re, im),
                Err(AnalyticCurveError::TauNotInUpperHalfPlane)
            );
            prop_assert!(!UpperHalfPlanePoint::is_in_upper_half_plane(&Complex64::new(re, im)));
        }
    }
}
