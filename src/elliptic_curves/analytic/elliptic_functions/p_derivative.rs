use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, ComplexLattice};

use crate::elliptic_curves::analytic::elliptic_functions::traits::{
    impl_elliptic_function_approximation, impl_has_pole_distance,
};
use crate::elliptic_curves::analytic::elliptic_functions::{
    EllipticFunctionTruncation, evaluator::evaluate_truncated_elliptic_function,
};

/// One truncated approximation to the derivative `℘′` of the Weierstrass
/// `℘`-function.
///
/// As with [`super::WeierstrassPApprox`], the stored value is a finite
/// educational approximation attached to a chosen truncation box in
/// lattice-index space.
#[derive(Clone, Debug, PartialEq)]
pub struct WeierstrassPDerivativeApprox {
    /// The original evaluation point supplied by the caller.
    z: Complex64,
    /// The truncated approximation to `℘′(z; Λ)`.
    value: Complex64,
    /// The truncation policy used for the lattice sum.
    truncation: EllipticFunctionTruncation,
    /// The number of nonzero lattice points used in the truncated sum.
    ///
    /// This count excludes the separate singular term `-2 / z³`.
    terms_used: usize,
    /// The smallest distance from the reduced evaluation point to the lattice
    /// poles inspected by this truncated approximation.
    pole_distance: f64,
}

impl WeierstrassPDerivativeApprox {
    pub(crate) fn from_parts(
        z: Complex64,
        value: Complex64,
        truncation: EllipticFunctionTruncation,
        terms_used: usize,
        pole_distance: f64,
    ) -> Self {
        Self {
            z,
            value,
            truncation,
            terms_used,
            pole_distance,
        }
    }

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

    /// Returns the smallest inspected distance to a lattice pole.
    pub fn pole_distance(&self) -> f64 {
        self.pole_distance
    }
}

impl_elliptic_function_approximation!(WeierstrassPDerivativeApprox);
impl_has_pole_distance!(WeierstrassPDerivativeApprox);

/// Approximates the derivative `℘′` of the Weierstrass `℘`-function attached
/// to `Λ`.
///
/// `℘′(z; Λ) = -2 / z³ - 2 Σ[1 / (z - ω)³]`
///
/// The infinite sum is truncated to the symmetric index box `-r ≤ m ≤ r`,
/// `-r ≤ n ≤ r`, omitting `(0, 0)`.
///
/// The input is first reduced modulo `Λ` to the chosen canonical representative.
/// Points numerically too close to a lattice pole are rejected with
/// [`AnalyticCurveError::PointTooCloseToLatticePoint`].
///
/// Complexity: `Θ(r²)` in the truncation radius `r`.
pub fn weierstrass_p_derivative(
    lattice: &ComplexLattice,
    z: Complex64,
    truncation: EllipticFunctionTruncation,
) -> Result<WeierstrassPDerivativeApprox, AnalyticCurveError> {
    evaluate_truncated_elliptic_function(
        lattice,
        z,
        truncation,
        |canonical_z| -2.0 * (Complex64::new(1.0, 0.0) / canonical_z.powu(3)),
        |canonical_z, omega| {
            let shifted = canonical_z - omega;

            -2.0 * (Complex64::new(1.0, 0.0) / shifted.powu(3))
        },
        WeierstrassPDerivativeApprox::from_parts,
    )
}
