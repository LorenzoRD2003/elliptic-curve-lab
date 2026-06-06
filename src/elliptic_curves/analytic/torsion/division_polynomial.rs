use crate::elliptic_curves::AffinePoint;
use crate::elliptic_curves::division_polynomials::{
    DivisionPolynomialXCriterionKind, division_polynomial_x_criterion_kind,
    evaluate_division_polynomial_x_criterion,
};
use crate::fields::ComplexApprox;

use crate::elliptic_curves::analytic::torsion::{
    AnalyticDivisionPolynomialComparisonCase, AnalyticDivisionPolynomialComparisonStatus,
    AnalyticEvenDivisionPolynomialReport, AnalyticOddDivisionPolynomialReport,
    EvenDivisionPolynomialVanishingBranch, map_torus_torsion_to_curve,
};
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice,
    EllipticFunctionTruncation, LatticeSumTruncation,
};

/// Compares analytic torus `n`-torsion against the complex division-polynomial
/// `x`-criterion on the associated analytic Weierstrass curve.
///
/// Each torus `n`-torsion class is first mapped to the curve via
/// `z ↦ (℘(z), ℘′(z))`, then the corresponding short-Weierstrass companion is
/// used to evaluate the complex division-polynomial `x`-criterion at
/// `x = ℘(z)`.
///
/// The identity class is reported separately as
/// [`AnalyticDivisionPolynomialComparisonStatus::PoleAtIdentity`], because it
/// maps to infinity and therefore has no finite `x`-coordinate.
///
/// The return type is case-split:
///
/// - [`AnalyticDivisionPolynomialComparisonCase::Pole`] for the identity class
/// - [`AnalyticDivisionPolynomialComparisonCase::Odd`] when the tested
///   criterion is `ψ_n(x)`
/// - [`AnalyticDivisionPolynomialComparisonCase::Even`] when the tested
///   criterion is `ε_n(x)` together with the active even-vanishing branch
///
/// Complexity: `Θ(n²(r_inv² + r_fun² + n^5))`, where the `n^5` term comes from
/// the current naive recursive division-polynomial construction for each
/// `x`-criterion evaluation.
pub fn compare_analytic_torsion_with_division_polynomial(
    lattice: &ComplexLattice,
    n: usize,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<Vec<AnalyticDivisionPolynomialComparisonCase>, AnalyticCurveError> {
    let criterion_kind =
        division_polynomial_x_criterion_kind(n).map_err(AnalyticCurveError::from)?;
    let curve = AnalyticWeierstrassCurve::from_lattice(lattice, invariant_truncation)?;
    let short_curve = curve.as_short_weierstrass();

    map_torus_torsion_to_curve(
        lattice,
        n,
        invariant_truncation,
        function_truncation,
        tolerance,
    )?
    .into_iter()
    .map(|torsion_point| match torsion_point.curve_point().clone() {
        AffinePoint::Infinity => Ok(AnalyticDivisionPolynomialComparisonCase::Pole {
            torsion_point,
            tolerance,
        }),
        AffinePoint::Finite { x, y } => {
            let value = evaluate_division_polynomial_x_criterion(&short_curve, n, &x)
                .map_err(AnalyticCurveError::from)?;
            let absolute = value.norm();
            let criterion_is_zero = ComplexApprox::is_zero_with_tolerance(&value, tolerance);

            match criterion_kind {
                DivisionPolynomialXCriterionKind::OddDivisionPolynomial => {
                    let status = if criterion_is_zero {
                        AnalyticDivisionPolynomialComparisonStatus::VanishesApproximately
                    } else {
                        AnalyticDivisionPolynomialComparisonStatus::DoesNotVanishApproximately
                    };

                    Ok(AnalyticDivisionPolynomialComparisonCase::Odd(
                        AnalyticOddDivisionPolynomialReport {
                            torsion_point,
                            x_value: x,
                            psi_n_x: value,
                            absolute_value: absolute,
                            status,
                            tolerance,
                        },
                    ))
                }
                DivisionPolynomialXCriterionKind::EvenYStrippedFactor => {
                    let y_is_zero = ComplexApprox::is_zero_with_tolerance(&y, tolerance);
                    let branch = match (y_is_zero, criterion_is_zero) {
                        (true, true) => EvenDivisionPolynomialVanishingBranch::BothBranches,
                        (true, false) => EvenDivisionPolynomialVanishingBranch::YApproxZero,
                        (false, true) => {
                            EvenDivisionPolynomialVanishingBranch::XCriterionApproxZero
                        }
                        (false, false) => EvenDivisionPolynomialVanishingBranch::NeitherBranch,
                    };
                    let status = if y_is_zero || criterion_is_zero {
                        AnalyticDivisionPolynomialComparisonStatus::VanishesApproximately
                    } else {
                        AnalyticDivisionPolynomialComparisonStatus::DoesNotVanishApproximately
                    };

                    Ok(AnalyticDivisionPolynomialComparisonCase::Even(
                        AnalyticEvenDivisionPolynomialReport {
                            torsion_point,
                            x_value: x,
                            epsilon_n_x: value,
                            absolute_value: absolute,
                            branch,
                            status,
                            tolerance,
                        },
                    ))
                }
            }
        }
    })
    .collect()
}

/// Restricts [`compare_analytic_torsion_with_division_polynomial`] to
/// primitive torus `n`-torsion classes.
///
/// Here “primitive” means exact torus order `n`, equivalently
/// `gcd(a, b, n) = 1`.
///
/// Complexity: `Θ(n²(r_inv² + r_fun² + n^5))`, since the current
/// implementation compares the full torus `n`-torsion grid first and then
/// filters.
pub fn compare_primitive_analytic_torsion_with_division_polynomial(
    lattice: &ComplexLattice,
    n: usize,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<Vec<AnalyticDivisionPolynomialComparisonCase>, AnalyticCurveError> {
    Ok(compare_analytic_torsion_with_division_polynomial(
        lattice,
        n,
        invariant_truncation,
        function_truncation,
        tolerance,
    )?
    .into_iter()
    .filter(|report| report.torsion_point().torus_point().index().is_primitive())
    .collect())
}
