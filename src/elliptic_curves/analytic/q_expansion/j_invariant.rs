use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ApproxTolerance, LatticeSumTruncation, UpperHalfPlanePoint,
    q_expansion::{
        JInvariantComparisonReport, ModularQExpansionCoefficients, ModularQParameter,
        QExpansionTruncation,
        family::{ModularQExpansionFamily, impl_modular_q_expansion_accessors},
    },
};
use crate::numerics::ComplexApproxComparison;

/// Marker type for the current educational truncated `q`-expansion family of
/// the modular `j`-invariant.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct JInvariantQExpansion;

/// Approximation report for the first truncated `q`-expansion of the modular
/// `j`-invariant.
///
/// Current convention: the stored truncation counts every retained coefficient
/// of the current educational table, starting at the principal term `q^{-1}`.
#[derive(Clone, Debug, PartialEq)]
pub struct JInvariantQExpansionApprox {
    q_parameter: ModularQParameter,
    value: Complex64,
    truncation: QExpansionTruncation,
}

impl JInvariantQExpansion {
    /// Returns the distinguished modular `j(q)` family object.
    pub fn new() -> Self {
        Self
    }

    /// Returns the currently shipped educational table for the modular
    /// `j`-invariant including the principal part:
    /// `q^{-1} + 744 + 196884q + 21493760q² + 864299970q³ + 20245856256q⁴`.
    pub(crate) fn current_table() -> ModularQExpansionCoefficients {
        ModularQExpansionCoefficients::from_integers(
            -1,
            vec![1, 744, 196_884, 21_493_760, 864_299_970, 20_245_856_256i64],
        )
    }

    /// Returns the number of coefficients currently shipped in the
    /// educational `j(q)` table, including the principal term `q^{-1}`.
    pub fn max_supported_terms() -> usize {
        Self::current_table().len()
    }

    /// Returns the richest truncation currently supported by the shipped
    /// educational `j(q)` table.
    pub fn full_current_table_truncation() -> QExpansionTruncation {
        QExpansionTruncation::new(Self::max_supported_terms())
            .expect("the shipped j(q) table always has positive length")
    }

    /// Returns the exact truncated coefficient table used for this `j(q)`
    /// approximation.
    #[allow(dead_code)]
    pub(crate) fn coefficients(
        truncation: QExpansionTruncation,
    ) -> Result<ModularQExpansionCoefficients, AnalyticCurveError> {
        let family = Self::new();
        <Self as ModularQExpansionFamily>::coefficients(&family, truncation)
    }

    /// Builds the truncated `j(q)` approximation from an upper-half-plane
    /// point `τ`.
    ///
    /// The current table is
    /// `j(q) = q^{-1} + 744 + 196884q + 21493760q² + 864299970q³ + 20245856256q⁴`.
    ///
    /// Because this is a `q`-expansion in the cusp coordinate, the
    /// approximation is especially effective when `|q|` is small, for example
    /// after reducing `τ` to a standard fundamental domain.
    pub fn evaluate_tau(
        tau: UpperHalfPlanePoint,
        truncation: QExpansionTruncation,
    ) -> Result<JInvariantQExpansionApprox, AnalyticCurveError> {
        let family = Self::new();
        <Self as ModularQExpansionFamily>::evaluate_tau(&family, tau, truncation)
    }

    /// Returns the currently shipped nonnegative-power coefficients of the
    /// modular `j`-invariant:
    /// `744, 196884, 21493760, 864299970, 20245856256`.
    ///
    /// The principal term `q^{-1}` is not included in this table.
    #[cfg(test)]
    pub(crate) fn non_negative_table() -> ModularQExpansionCoefficients {
        ModularQExpansionCoefficients::from_integers(
            0,
            vec![744, 196_884, 21_493_760, 864_299_970, 20_245_856_256i64],
        )
    }

    /// Compares the two current analytic routes to the modular `j`-invariant:
    ///
    /// - truncated Eisenstein sums on the lattice `Λ_τ = ℤ + ℤτ`
    /// - truncated cusp expansion in `q = e^{2π i τ}`
    ///
    /// This is an educational numerical experiment rather than a certified
    /// modular-forms routine. Its quality depends both on the lattice truncation
    /// radius and on how small `|q|` is for the chosen `τ`.
    ///
    /// Complexity: `Θ(r² + N)`, where `r` is the lattice truncation radius and
    /// `N = q_truncation.terms()`.
    pub fn compare_with_eisenstein_sum(
        tau: UpperHalfPlanePoint,
        lattice_truncation: LatticeSumTruncation,
        q_truncation: QExpansionTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<JInvariantComparisonReport, AnalyticCurveError> {
        let invariants = tau.analytic_invariants(lattice_truncation)?;
        let q_approximation = JInvariantQExpansion::evaluate_tau(tau.clone(), q_truncation)?;
        Ok(JInvariantComparisonReport::new(
            tau,
            ComplexApproxComparison::new(
                *invariants.j_invariant(),
                *q_approximation.value(),
                tolerance,
            ),
            lattice_truncation,
            q_truncation,
        ))
    }
}

impl ModularQExpansionFamily for JInvariantQExpansion {
    type Approximation = JInvariantQExpansionApprox;

    fn coefficients(
        &self,
        truncation: QExpansionTruncation,
    ) -> Result<ModularQExpansionCoefficients, AnalyticCurveError> {
        Self::current_table().truncated(truncation)
    }

    fn build_approximation(
        &self,
        q_parameter: ModularQParameter,
        value: Complex64,
        truncation: QExpansionTruncation,
    ) -> Self::Approximation {
        JInvariantQExpansionApprox {
            q_parameter,
            value,
            truncation,
        }
    }
}

impl_modular_q_expansion_accessors!(JInvariantQExpansionApprox);
