use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticCurveMembershipReport, AnalyticCurvePoint,
    AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice, EllipticFunctionTruncation,
    FundamentalParallelogramCoordinate, LatticeSumTruncation, weierstrass_p,
    weierstrass_p_derivative,
};
use crate::fields::ComplexApprox;

/// Finite or pole-valued output of the torus-to-curve map `z ↦ (℘(z), ℘′(z))`.
///
/// At lattice points, the classical meromorphic coordinates `℘(z)` and
/// `℘′(z)` blow up. Geometrically, those torus points map to the point at
/// infinity on the corresponding analytic Weierstrass curve, so this enum
/// records either that pole case or the finite pair of evaluated coordinates.
#[derive(Clone, Debug, PartialEq)]
pub enum TorusToCurveValues {
    /// The input represents a lattice point, so the map lands at infinity.
    Pole,
    /// The input has finite `℘` and `℘′` values.
    FiniteValues { p: Complex64, p_prime: Complex64 },
}

/// Structured result of mapping one torus point to the analytic Weierstrass
/// curve attached to the same lattice.
#[derive(Clone, Debug, PartialEq)]
pub struct TorusToCurveMapResult {
    z: Complex64,
    point: AnalyticCurvePoint,
    curve: AnalyticWeierstrassCurve,
    values: TorusToCurveValues,
    membership_report: AnalyticCurveMembershipReport,
}

impl TorusToCurveMapResult {
    /// Returns the original complex torus representative requested by the caller.
    pub fn z(&self) -> &Complex64 {
        &self.z
    }

    /// Returns the resulting analytic curve point.
    pub fn point(&self) -> &AnalyticCurvePoint {
        &self.point
    }

    /// Returns the analytic Weierstrass curve used for the map.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the finite `℘`/`℘′` pair or the pole case.
    pub fn values(&self) -> &TorusToCurveValues {
        &self.values
    }

    /// Returns the approximate curve-membership report for `point`.
    pub fn membership_report(&self) -> &AnalyticCurveMembershipReport {
        &self.membership_report
    }

    /// Returns whether the mapped point was accepted as lying on the curve
    /// under the caller-provided tolerance.
    pub fn lies_on_curve(&self) -> bool {
        self.membership_report.is_on_curve()
    }
}

/// Maps one complex torus point to the analytic Weierstrass curve attached to
/// the same lattice.
///
/// This implements the analytic correspondence `z ↦ (℘(z), ℘′(z))`
/// between the complex torus `ℂ / Λ` and the cubic `y² = 4x³ - g₂x - g₃`.
///
/// If `z` represents a lattice point, then `℘(z)` and `℘′(z)` are poles and
/// the map lands at the point at infinity instead of returning an error.
///
/// Complexity: `Θ(r_inv² + r_fun²)`, where `r_inv` is the invariant
/// truncation radius and `r_fun` is the elliptic-function truncation radius.
pub fn map_torus_point_to_curve(
    lattice: &ComplexLattice,
    z: Complex64,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<TorusToCurveMapResult, AnalyticCurveError> {
    let reduced_coordinate = lattice.reduce_complex_point_to_fundamental_coordinates(z)?;
    let reduced_z = lattice.point_from_fundamental_coordinates(reduced_coordinate);
    let curve = AnalyticWeierstrassCurve::from_lattice(lattice, invariant_truncation)?;

    // reduced_z == 0 is equivalent to z being a lattice point
    if ComplexApprox::is_zero_with_tolerance(&reduced_z, ComplexApprox::default_tolerance()) {
        let point = AnalyticCurvePoint::infinity();
        let membership_report = curve.membership_report(&point, tolerance);

        return Ok(TorusToCurveMapResult {
            z,
            point,
            curve,
            values: TorusToCurveValues::Pole,
            membership_report,
        });
    }

    let p_approx = weierstrass_p(lattice, z, function_truncation)?;
    let p_prime_approx = weierstrass_p_derivative(lattice, z, function_truncation)?;
    let point = AnalyticCurvePoint::new(*p_approx.value(), *p_prime_approx.value());
    let membership_report = curve.membership_report(&point, tolerance);

    Ok(TorusToCurveMapResult {
        z,
        point: point.clone(),
        curve,
        values: TorusToCurveValues::FiniteValues {
            p: *p_approx.value(),
            p_prime: *p_prime_approx.value(),
        },
        membership_report,
    })
}

/// Maps one canonical fundamental-parallelogram coordinate to the analytic
/// Weierstrass curve attached to the same lattice.
///
/// Complexity: `Θ(r_inv² + r_fun²)`, inherited from [`map_torus_point_to_curve`].
pub fn map_fundamental_point_to_curve(
    lattice: &ComplexLattice,
    coord: FundamentalParallelogramCoordinate,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<TorusToCurveMapResult, AnalyticCurveError> {
    let z = lattice.point_from_fundamental_coordinates(coord);
    map_torus_point_to_curve(
        lattice,
        z,
        invariant_truncation,
        function_truncation,
        tolerance,
    )
}
