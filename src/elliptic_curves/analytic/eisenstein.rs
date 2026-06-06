use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, ComplexLattice, LatticeSumTruncation};
use crate::numerics::ComplexDifferenceReport;

/// Finite square-box approximation to the lattice Eisenstein sum `G_k(Λ)`.
///
/// The reported value approximates
/// `G_k(Λ) = Σ_{ω ∈ Λ \ {0}} 1 / ω^k`
/// by restricting to the punctured index box
/// `-r ≤ m ≤ r`, `-r ≤ n ≤ r`, `(m, n) ≠ (0, 0)`.
///
/// The approximation keeps the truncation policy and the number of nonzero
/// lattice terms explicit so later educational reports can compare truncation
/// scales directly.
#[derive(Clone, Debug, PartialEq)]
pub struct EisensteinSumApprox {
    /// Weight `k` in `G_k(Λ)`.
    pub weight: usize,
    /// Truncated lattice-sum value in `ℂ`.
    pub value: Complex64,
    /// Square-box truncation used for the approximation.
    pub truncation: LatticeSumTruncation,
    /// Number of nonzero lattice terms conceptually included in the box.
    pub terms_used: usize,
}

/// Side-by-side comparison between two square-box truncations of the same
/// lattice Eisenstein sum.
///
/// This report does not certify convergence in a theorem-proof sense.
#[derive(Clone, Debug, PartialEq)]
pub struct TruncationConvergenceReport {
    /// Smaller square-box truncation used in the comparison.
    small: LatticeSumTruncation,
    /// Larger square-box truncation used in the comparison.
    large: LatticeSumTruncation,
    /// Shared left/right residual payload for the two truncations.
    comparison: ComplexDifferenceReport,
}

impl TruncationConvergenceReport {
    /// Returns the smaller truncation radius.
    pub fn small(&self) -> LatticeSumTruncation {
        self.small
    }

    /// Returns the larger truncation radius.
    pub fn large(&self) -> LatticeSumTruncation {
        self.large
    }

    /// Returns the truncated value computed with `small`.
    pub fn small_value(&self) -> &Complex64 {
        self.comparison.right()
    }

    /// Returns the truncated value computed with `large`.
    pub fn large_value(&self) -> &Complex64 {
        self.comparison.left()
    }

    /// Returns the residual `large_value - small_value`.
    pub fn difference(&self) -> &Complex64 {
        self.comparison.difference()
    }

    /// Returns the Euclidean norm `|large_value - small_value|`.
    pub fn absolute_difference(&self) -> f64 {
        self.comparison.absolute_difference()
    }
}

/// Approximates the lattice Eisenstein sum `G_k(Λ)` by a finite punctured
/// square-box sum. It uses the truncation
/// `-r ≤ m ≤ r`, `-r ≤ n ≤ r`, `(m, n) ≠ (0, 0)` and evaluates
/// `Σ 1 / (mω₁ + nω₂)^k` directly.
///
/// Only weights `k ≥ 3` are accepted. That threshold matches the
/// absolute-convergence condition for the raw lattice sum in two dimensions.
/// The low weights `k = 0, 1, 2` are rejected because the naive lattice sum
/// does not define an absolutely convergent Eisenstein series there.
///
/// For odd weights `k > 2`, the function returns exactly `0` without
/// enumerating the lattice box. This uses the symmetry `ω ↔ -ω` in `Λ`.
///
/// Complexity: for even weights, `Θ(r²)` time in the truncation radius `r`
/// and `Θ(r²)` memory through the current boxed lattice-point enumeration.
/// For odd weights `k > 2`, the current implementation short-circuits in
/// `Θ(1)`.
pub fn eisenstein_sum(
    lattice: &ComplexLattice,
    weight: usize,
    truncation: LatticeSumTruncation,
) -> Result<EisensteinSumApprox, AnalyticCurveError> {
    if weight < 3 {
        return Err(AnalyticCurveError::InvalidEisensteinWeight);
    }

    let terms_used = truncation.nonzero_terms_in_square_box();

    if weight % 2 == 1 {
        return Ok(EisensteinSumApprox {
            weight,
            value: Complex64::new(0.0, 0.0),
            truncation,
            terms_used,
        });
    }

    let value = lattice
        .nonzero_lattice_points_in_box(truncation.radius())
        .into_iter()
        .fold(Complex64::new(0.0, 0.0), |acc, point| {
            acc + Complex64::new(1.0, 0.0) / point.value.powu(weight as u32)
        });

    Ok(EisensteinSumApprox {
        weight,
        value,
        truncation,
        terms_used,
    })
}

/// Approximates the classical weight-4 Eisenstein sum `G₄(Λ)`.
///
/// Complexity: `Θ(r²)` in the truncation radius `r`.
pub fn g4_sum(
    lattice: &ComplexLattice,
    truncation: LatticeSumTruncation,
) -> Result<EisensteinSumApprox, AnalyticCurveError> {
    eisenstein_sum(lattice, 4, truncation)
}

/// Approximates the classical weight-6 Eisenstein sum `G₆(Λ)`.
///
/// Complexity: `Θ(r²)` in the truncation radius `r`.
pub fn g6_sum(
    lattice: &ComplexLattice,
    truncation: LatticeSumTruncation,
) -> Result<EisensteinSumApprox, AnalyticCurveError> {
    eisenstein_sum(lattice, 6, truncation)
}

/// Compares two truncations of the same lattice Eisenstein sum.
///
/// The argument named `large` must have strictly larger radius than `small`.
/// If not, the comparison request is rejected as malformed.
///
/// The returned report records both truncated values together with the
/// difference `large_value - small_value` and its absolute magnitude.
///
/// Complexity: `Θ(r_small² + r_large²)` time, since both truncations are
/// evaluated independently.
pub fn compare_eisenstein_truncations(
    lattice: &ComplexLattice,
    weight: usize,
    small: LatticeSumTruncation,
    large: LatticeSumTruncation,
) -> Result<TruncationConvergenceReport, AnalyticCurveError> {
    if large.radius() <= small.radius() {
        return Err(AnalyticCurveError::InvalidTruncationComparison);
    }

    let small_sum = eisenstein_sum(lattice, weight, small)?;
    let large_sum = eisenstein_sum(lattice, weight, large)?;
    Ok(TruncationConvergenceReport {
        small,
        large,
        comparison: ComplexDifferenceReport::new(large_sum.value, small_sum.value),
    })
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::elliptic_curves::analytic::{
        AnalyticCurveError, ComplexLattice, LatticeSumTruncation, UpperHalfPlanePoint,
    };
    use crate::elliptic_curves::analytic::{
        EisensteinSumApprox, compare_eisenstein_truncations, eisenstein_sum, g4_sum, g6_sum,
    };
    use crate::fields::{ComplexApprox, Field};

    fn standard_square_lattice() -> ComplexLattice {
        ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
    }

    fn standard_hexagonal_lattice() -> ComplexLattice {
        ComplexLattice::from_tau(UpperHalfPlanePoint::tau_rho())
    }

    #[test]
    fn eisenstein_sum_rejects_non_convergent_low_weights() {
        let lattice = standard_square_lattice();
        let truncation = LatticeSumTruncation::default_educational();

        assert_eq!(
            eisenstein_sum(&lattice, 0, truncation),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
        assert_eq!(
            eisenstein_sum(&lattice, 1, truncation),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
        assert_eq!(
            eisenstein_sum(&lattice, 2, truncation),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
    }

    #[test]
    fn odd_weight_eisenstein_sum_short_circuits_to_exact_zero() {
        let lattice = standard_hexagonal_lattice();
        let truncation = LatticeSumTruncation::larger_for_comparison();

        let sum = eisenstein_sum(&lattice, 5, truncation).unwrap();

        assert_eq!(
            sum,
            EisensteinSumApprox {
                weight: 5,
                value: Complex64::new(0.0, 0.0),
                truncation,
                terms_used: truncation.nonzero_terms_in_square_box(),
            }
        );
    }

    #[test]
    fn g4_sum_matches_direct_even_weight_dispatch() {
        let lattice = standard_square_lattice();
        let truncation = LatticeSumTruncation::default_educational();

        let generic = eisenstein_sum(&lattice, 4, truncation).unwrap();
        let specialized = g4_sum(&lattice, truncation).unwrap();

        assert_eq!(generic, specialized);
        assert_eq!(generic.weight, 4);
        assert_eq!(generic.terms_used, truncation.nonzero_terms_in_square_box());
    }

    #[test]
    fn g6_sum_matches_direct_even_weight_dispatch() {
        let lattice = standard_hexagonal_lattice();
        let truncation = LatticeSumTruncation::default_educational();

        let generic = eisenstein_sum(&lattice, 6, truncation).unwrap();
        let specialized = g6_sum(&lattice, truncation).unwrap();

        assert_eq!(generic, specialized);
        assert_eq!(generic.weight, 6);
        assert_eq!(generic.terms_used, truncation.nonzero_terms_in_square_box());
    }

    #[test]
    fn square_lattice_g4_truncation_is_real_and_positive() {
        let lattice = standard_square_lattice();
        let truncation = LatticeSumTruncation::default_educational();

        let sum = g4_sum(&lattice, truncation).unwrap();

        assert!(sum.value.re > 0.0);
        assert!(ComplexApprox::default_tolerance().real_close(sum.value.im, 0.0));
    }

    #[test]
    fn larger_even_weight_truncation_changes_the_approximation() {
        let lattice = standard_hexagonal_lattice();
        let small = g6_sum(&lattice, LatticeSumTruncation::default_educational()).unwrap();
        let larger = g6_sum(&lattice, LatticeSumTruncation::larger_for_comparison()).unwrap();

        assert_eq!(small.weight, 6);
        assert_eq!(larger.weight, 6);
        assert!(larger.terms_used > small.terms_used);
        assert!(!ComplexApprox::eq(&small.value, &larger.value));
    }

    #[test]
    fn truncation_comparison_requires_strictly_increasing_radius() {
        let lattice = standard_square_lattice();
        let same = LatticeSumTruncation::default_educational();
        let larger = LatticeSumTruncation::larger_for_comparison();

        assert_eq!(
            compare_eisenstein_truncations(&lattice, 4, same, same),
            Err(AnalyticCurveError::InvalidTruncationComparison)
        );
        assert_eq!(
            compare_eisenstein_truncations(&lattice, 4, larger, same),
            Err(AnalyticCurveError::InvalidTruncationComparison)
        );
    }

    #[test]
    fn truncation_comparison_reports_both_values_and_difference() {
        let lattice = standard_square_lattice();
        let small = LatticeSumTruncation::default_educational();
        let large = LatticeSumTruncation::larger_for_comparison();

        let report = compare_eisenstein_truncations(&lattice, 4, small, large).unwrap();
        let small_sum = eisenstein_sum(&lattice, 4, small).unwrap();
        let large_sum = eisenstein_sum(&lattice, 4, large).unwrap();
        let expected_difference = large_sum.value - small_sum.value;

        assert_eq!(report.small(), small);
        assert_eq!(report.large(), large);
        assert_eq!(report.small_value(), &small_sum.value);
        assert_eq!(report.large_value(), &large_sum.value);
        assert_eq!(report.difference(), &expected_difference);
        assert_eq!(report.absolute_difference(), expected_difference.norm());
    }

    #[test]
    fn odd_weight_truncation_comparison_is_exactly_zero_on_both_sides() {
        let lattice = standard_hexagonal_lattice();
        let small = LatticeSumTruncation::default_educational();
        let large = LatticeSumTruncation::larger_for_comparison();

        let report = compare_eisenstein_truncations(&lattice, 5, small, large).unwrap();

        assert_eq!(report.small_value(), &Complex64::new(0.0, 0.0));
        assert_eq!(report.large_value(), &Complex64::new(0.0, 0.0));
        assert_eq!(report.difference(), &Complex64::new(0.0, 0.0));
        assert_eq!(report.absolute_difference(), 0.0);
    }
}
