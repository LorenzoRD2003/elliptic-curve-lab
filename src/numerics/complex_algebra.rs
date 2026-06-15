use num_complex::Complex64;

use crate::numerics::ApproxTolerance;

/// Returns the principal complex cube root of `z`.
pub(crate) fn principal_cube_root(z: Complex64) -> Complex64 {
    let radius = z.norm().cbrt();
    let angle = z.arg() / 3.0;
    Complex64::from_polar(radius, angle)
}

/// Returns a fixed primitive cube root of unity
/// `ω = -1/2 + (√3/2)i`.
pub(crate) fn primitive_cube_root_of_unity() -> Complex64 {
    Complex64::new(-0.5, f64::sqrt(3.0) * 0.5)
}

/// Returns the three complex cube-root branches of `z`, in principal-branch
/// order `{ζ, ωζ, ω²ζ}`.
pub(crate) fn cube_root_branches(z: Complex64) -> [Complex64; 3] {
    if z == Complex64::new(0.0, 0.0) {
        return [Complex64::new(0.0, 0.0); 3];
    }

    let principal = principal_cube_root(z);
    let omega = primitive_cube_root_of_unity();
    let omega_sq = omega * omega;

    [principal, omega * principal, omega_sq * principal]
}

/// Heuristically detects when the depressed cubic `x^3 + px + q = 0`
/// is numerically close to its pure-cubic limit `x^3 + q = 0`.
pub(crate) fn is_near_pure_cubic_regime(
    p: Complex64,
    q: Complex64,
    tolerance: ApproxTolerance,
) -> bool {
    let tolerance_scale = tolerance.absolute.max(tolerance.relative).max(f64::EPSILON);
    let q_scale_squared = q.norm().cbrt().powi(2).max(tolerance_scale);
    let normalized_p = p.norm() / q_scale_squared;
    normalized_p <= tolerance_scale.powf(0.25)
}
