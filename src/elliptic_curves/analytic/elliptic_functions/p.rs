use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, ComplexLattice};

use super::traits::{impl_elliptic_function_approximation, impl_has_pole_distance};
use super::{EllipticFunctionTruncation, evaluator::evaluate_truncated_elliptic_function};

/// One truncated approximation to the Weierstrass `℘`-function.
///
/// The stored value is an educational finite approximation to `℘(z; Λ)`
/// obtained from a symmetric square-box lattice truncation. The input `z` is
/// stored exactly as it was requested, even though the implementation may
/// first reduce it modulo `Λ` internally to choose a canonical representative
/// in the fundamental parallelogram.
#[derive(Clone, Debug, PartialEq)]
pub struct WeierstrassPApprox {
    /// The original evaluation point supplied by the caller.
    z: Complex64,
    /// The truncated approximation to `℘(z; Λ)`.
    value: Complex64,
    /// The truncation policy used for the lattice sum.
    truncation: EllipticFunctionTruncation,
    /// The number of nonzero lattice points used in the truncated sum.
    ///
    /// This count excludes the separate singular term `1 / z²`.
    terms_used: usize,
    /// The smallest distance from the reduced evaluation point to the lattice
    /// poles inspected by this truncated approximation.
    pole_distance: f64,
}

impl_elliptic_function_approximation!(WeierstrassPApprox);
impl_has_pole_distance!(WeierstrassPApprox);

/// Approximates the Weierstrass `℘`-function attached to `Λ`.
///
/// `℘(z; Λ) = 1 / z² + Σ_{ω ∈ Λ \ {0}} [(1 / (z - ω)²) - (1 / ω²)]`
///
/// The infinite sum is truncated to the symmetric index box `-r ≤ m ≤ r`,
/// `-r ≤ n ≤ r`, omitting `(0, 0)`.
///
/// Because `℘` is periodic modulo `Λ`, the computation first reduces `z` to
/// the chosen half-open fundamental parallelogram. If that representative is
/// numerically too close to `0`, the evaluation is rejected with
/// [`AnalyticCurveError::PointTooCloseToLatticePoint`], since `℘` has poles
/// exactly at lattice points.
///
/// Complexity: `Θ(r²)` in the truncation radius `r`.
pub fn weierstrass_p(
    lattice: &ComplexLattice,
    z: Complex64,
    truncation: EllipticFunctionTruncation,
) -> Result<WeierstrassPApprox, AnalyticCurveError> {
    evaluate_truncated_elliptic_function(
        lattice,
        z,
        truncation,
        |canonical_z| Complex64::new(1.0, 0.0) / canonical_z.powu(2),
        |canonical_z, omega| {
            let one = Complex64::new(1.0, 0.0);
            let shifted = canonical_z - omega;

            (one / shifted.powu(2)) - (one / omega.powu(2))
        },
        |z, value, truncation, terms_used, pole_distance| WeierstrassPApprox {
            z,
            value,
            truncation,
            terms_used,
            pole_distance,
        },
    )
}
