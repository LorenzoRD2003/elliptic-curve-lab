use num_complex::Complex64;

use crate::fields::ComplexApprox;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, EllipticFunctionTruncation, HasPoleDistance,
};

/// One truncated approximation to the Weierstrass `ζ`-function.
///
/// Unlike `℘` and `℘′`, the Weierstrass `ζ`-function is not elliptic: it is
/// only quasi-periodic. So the stored value approximates the actual meromorphic
/// function on `ℂ`, not a torus-descended function on `ℂ / Λ`.
#[derive(Clone, Debug, PartialEq)]
pub struct WeierstrassZetaApprox {
    /// The original evaluation point supplied by the caller.
    z: Complex64,
    /// The truncated approximation to `ζ(z; Λ)`.
    value: Complex64,
    /// The truncation policy used for the lattice sum.
    truncation: EllipticFunctionTruncation,
    /// The number of nonzero lattice points used in the truncated sum.
    ///
    /// This count excludes the separate singular term `1 / z`.
    terms_used: usize,
    /// The smallest distance from the reduced evaluation point to the lattice
    /// poles inspected by this truncated approximation.
    pole_distance: f64,
}

impl WeierstrassZetaApprox {
    /// Returns the original evaluation point supplied by the caller.
    pub fn z(&self) -> &Complex64 {
        &self.z
    }

    /// Returns the approximate complex value produced by the truncation.
    pub fn value(&self) -> &Complex64 {
        &self.value
    }

    /// Returns the truncation policy used for this approximation.
    pub fn truncation(&self) -> EllipticFunctionTruncation {
        self.truncation
    }

    /// Returns the number of nonzero lattice terms that were summed.
    pub fn terms_used(&self) -> usize {
        self.terms_used
    }
}

impl HasPoleDistance for WeierstrassZetaApprox {
    fn pole_distance(&self) -> f64 {
        self.pole_distance
    }
}

/// Approximates the Weierstrass `ζ`-function attached to `Λ`.
///
/// `ζ(z; Λ) = 1 / z + Σ_{ω ∈ Λ \ {0}} [1 / (z - ω) + 1 / ω + z / ω²]`
///
/// The infinite sum is truncated to the symmetric index box `-r ≤ m ≤ r`,
/// `-r ≤ n ≤ r`, omitting `(0, 0)`.
///
/// A key subtlety is that `ζ` is not elliptic, only quasi-periodic. So unlike
/// `℘` and `℘′`, the evaluation must use the original point `z` itself rather
/// than first reducing modulo `Λ`. We still use the reduced representative only
/// to detect whether `z` lies numerically too close to some lattice pole.
///
/// Complexity: `Θ(r²)` time and `Θ(r²)` memory in the truncation radius `r`,
/// because the current implementation enumerates the full nonzero square box
/// before summing.
pub fn weierstrass_zeta(
    lattice: &ComplexLattice,
    z: Complex64,
    truncation: EllipticFunctionTruncation,
) -> Result<WeierstrassZetaApprox, AnalyticCurveError> {
    let canonical_coordinate = lattice.reduce_complex_point_to_fundamental_coordinates(z)?;
    let canonical_z = lattice.point_from_fundamental_coordinates(canonical_coordinate);
    let tolerance = ComplexApprox::default_tolerance();

    if ComplexApprox::is_zero_with_tolerance(&canonical_z, tolerance) {
        return Err(AnalyticCurveError::PointTooCloseToLatticePoint);
    }

    let lattice_points = lattice.nonzero_lattice_points_in_box(truncation.radius());
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

    Ok(WeierstrassZetaApprox {
        z,
        value,
        truncation,
        terms_used: lattice_points.len(),
        pole_distance,
    })
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::elliptic_curves::analytic::{WeierstrassZetaApprox, weierstrass_zeta};
    use crate::{
        elliptic_curves::analytic::{
            AnalyticCurveError, ComplexLattice, EllipticFunctionApproximation,
            EllipticFunctionTruncation, HasPoleDistance, UpperHalfPlanePoint, weierstrass_p,
        },
        fields::ComplexApprox,
    };

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    fn square_lattice() -> ComplexLattice {
        ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
    }

    #[test]
    fn weierstrass_zeta_rejects_lattice_points_as_poles() {
        let lattice = square_lattice();
        let truncation = EllipticFunctionTruncation::default_educational();

        assert_eq!(
            weierstrass_zeta(&lattice, c(0.0, 0.0), truncation),
            Err(AnalyticCurveError::PointTooCloseToLatticePoint)
        );
        assert_eq!(
            weierstrass_zeta(&lattice, c(1.0, 0.0), truncation),
            Err(AnalyticCurveError::PointTooCloseToLatticePoint)
        );
    }

    #[test]
    fn weierstrass_zeta_reports_input_truncation_and_terms_used() {
        let lattice = square_lattice();
        let truncation = EllipticFunctionTruncation::default_educational();
        let z = c(0.3, 0.2);

        let approximation = weierstrass_zeta(&lattice, z, truncation).unwrap();

        assert_eq!(
            approximation,
            WeierstrassZetaApprox {
                z,
                value: *approximation.value(),
                truncation,
                terms_used: 24,
                pole_distance: approximation.pole_distance(),
            }
        );
    }

    #[test]
    fn weierstrass_zeta_is_odd() {
        let lattice = square_lattice();
        let truncation = EllipticFunctionTruncation::default_educational();
        let z = c(0.31, 0.22);

        let positive = weierstrass_zeta(&lattice, z, truncation).unwrap();
        let negative = weierstrass_zeta(&lattice, -z, truncation).unwrap();

        assert!(ComplexApprox::eq_with_tolerance(
            negative.value(),
            &(-*positive.value()),
            ComplexApprox::default_tolerance()
        ));
    }

    #[test]
    fn weierstrass_zeta_derivative_matches_minus_weierstrass_p_by_finite_difference() {
        let lattice = square_lattice();
        let truncation = EllipticFunctionTruncation::default_educational();
        let z = c(0.23, 0.19);
        let h = 1.0e-6;

        let p = weierstrass_p(&lattice, z, truncation).unwrap();
        let zeta_forward = weierstrass_zeta(&lattice, z + c(h, 0.0), truncation).unwrap();
        let zeta_backward = weierstrass_zeta(&lattice, z - c(h, 0.0), truncation).unwrap();
        let finite_difference = (*zeta_forward.value() - *zeta_backward.value()) / (2.0 * h);

        assert!(ComplexApprox::eq_with_tolerance(
            &finite_difference,
            &(-*p.value()),
            crate::numerics::ApproxTolerance::new(1.0e-5, 1.0e-5)
        ));
    }
}
