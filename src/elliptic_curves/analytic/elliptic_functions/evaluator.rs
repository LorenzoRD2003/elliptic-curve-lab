use num_complex::Complex64;

use crate::{
    elliptic_curves::analytic::{AnalyticCurveError, ComplexLattice},
    fields::ComplexApprox,
};

use super::EllipticFunctionTruncation;

/// Shared reduction-and-summation routine for truncated elliptic-function
/// lattice expressions.
///
/// This helper is meant for callers who want to build a first custom elliptic
/// function from the same basic pattern used by `℘` and `℘′`:
/// - reduce `z` modulo `Λ`
/// - reject poles at lattice points
/// - add one singular term at the reduced representative
/// - add one lattice term for each nonzero `ω` in the chosen square box
/// - package the result into any caller-chosen report type
///
/// Complexity: `Θ(r²)` time and `Θ(r²)` memory in the truncation radius `r`,
/// because the current implementation enumerates the full nonzero square box
/// before summing.
pub fn evaluate_truncated_elliptic_function<A, S, L, B>(
    lattice: &ComplexLattice,
    z: Complex64,
    truncation: EllipticFunctionTruncation,
    singular_term: S,
    lattice_term: L,
    build_approximation: B,
) -> Result<A, AnalyticCurveError>
where
    S: Fn(Complex64) -> Complex64,
    L: Fn(Complex64, Complex64) -> Complex64,
    B: FnOnce(Complex64, Complex64, EllipticFunctionTruncation, usize, f64) -> A,
{
    let canonical_coordinate = lattice.reduce_complex_point_to_fundamental_coordinates(z)?;
    let canonical_z = lattice.point_from_fundamental_coordinates(canonical_coordinate);
    let tolerance = ComplexApprox::default_tolerance();

    if ComplexApprox::is_zero_with_tolerance(&canonical_z, tolerance) {
        return Err(AnalyticCurveError::PointTooCloseToLatticePoint);
    }

    let lattice_points = lattice.nonzero_lattice_points_in_box(truncation.radius());
    let mut value = singular_term(canonical_z);
    let mut pole_distance = canonical_z.norm();

    for point in &lattice_points {
        let shifted = canonical_z - point.value;
        let shifted_norm = shifted.norm();

        if ComplexApprox::is_zero_with_tolerance(&shifted, tolerance) {
            return Err(AnalyticCurveError::PointTooCloseToLatticePoint);
        }

        pole_distance = pole_distance.min(shifted_norm);
        value += lattice_term(canonical_z, point.value);
    }

    if !value.re.is_finite() || !value.im.is_finite() {
        return Err(AnalyticCurveError::NumericalComparisonFailed);
    }

    Ok(build_approximation(
        z,
        value,
        truncation,
        lattice_points.len(),
        pole_distance,
    ))
}
