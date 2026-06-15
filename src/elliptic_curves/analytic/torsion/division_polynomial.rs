use crate::elliptic_curves::{
    AffinePoint,
    analytic::{
        AnalyticCurveError, AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice,
        EllipticFunctionTruncation, LatticeSumTruncation,
        torsion::{
            AnalyticDivisionPolynomialComparisonCase, AnalyticDivisionPolynomialComparisonStatus,
            AnalyticEvenDivisionPolynomialReport, AnalyticOddDivisionPolynomialReport,
            EvenDivisionPolynomialVanishingBranch,
        },
    },
};
use crate::fields::complex_approx::ComplexApprox;

impl ComplexLattice {
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
        &self,
        n: usize,
        invariant_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<Vec<AnalyticDivisionPolynomialComparisonCase>, AnalyticCurveError> {
        let curve = AnalyticWeierstrassCurve::from_lattice(self, invariant_truncation)?;
        let short_curve = curve.as_short_weierstrass();
        let uses_y_factor = short_curve
            .division_polynomial_uses_y_factor(n)
            .map_err(AnalyticCurveError::from)?;

        self.map_torus_torsion_to_curve(n, invariant_truncation, function_truncation, tolerance)?
            .into_iter()
            .map(|torsion_point| match torsion_point.curve_point().clone() {
                AffinePoint::Infinity => Ok(AnalyticDivisionPolynomialComparisonCase::Pole {
                    torsion_point,
                    tolerance,
                }),
                AffinePoint::Finite { x, y } => {
                    let value = short_curve
                        .evaluate_division_polynomial_x_criterion(n, &x)
                        .map_err(AnalyticCurveError::from)?;
                    let absolute = value.norm();
                    let criterion_is_zero =
                        ComplexApprox::is_zero_with_tolerance(&value, tolerance);

                    if !uses_y_factor {
                        let status = if criterion_is_zero {
                            AnalyticDivisionPolynomialComparisonStatus::VanishesApproximately
                        } else {
                            AnalyticDivisionPolynomialComparisonStatus::DoesNotVanishApproximately
                        };

                        Ok(AnalyticDivisionPolynomialComparisonCase::Odd(
                            AnalyticOddDivisionPolynomialReport::new(
                                torsion_point,
                                x,
                                value,
                                absolute,
                                status,
                                tolerance,
                            ),
                        ))
                    } else {
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
                            AnalyticEvenDivisionPolynomialReport::new(
                                torsion_point,
                                x,
                                value,
                                absolute,
                                branch,
                                status,
                                tolerance,
                            ),
                        ))
                    }
                }
            })
            .collect()
    }

    /// Restricts [`Self::compare_analytic_torsion_with_division_polynomial`] to
    /// primitive torus `n`-torsion classes.
    ///
    /// Here “primitive” means exact torus order `n`, equivalently
    /// `gcd(a, b, n) = 1`.
    ///
    /// Complexity: `Θ(n²(r_inv² + r_fun² + n^5))`, since the current
    /// implementation compares the full torus `n`-torsion grid first and then
    /// filters.
    pub fn compare_primitive_analytic_torsion_with_division_polynomial(
        &self,
        n: usize,
        invariant_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<Vec<AnalyticDivisionPolynomialComparisonCase>, AnalyticCurveError> {
        Ok(self
            .compare_analytic_torsion_with_division_polynomial(
                n,
                invariant_truncation,
                function_truncation,
                tolerance,
            )?
            .into_iter()
            .filter(|report| report.torsion_point().index().is_primitive())
            .collect())
    }
}
