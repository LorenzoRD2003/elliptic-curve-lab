use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, UpperHalfPlanePoint};

use crate::elliptic_curves::analytic::q_expansion::family::{
    ModularQExpansionFamily, impl_modular_q_expansion_accessors,
};
use crate::elliptic_curves::analytic::q_expansion::{
    ModularQExpansionCoefficients, ModularQParameter, QExpansionTruncation,
};

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

    /// Returns the full currently shipped coefficient table for the
    /// educational `j(q)` experiment.
    pub fn current_table() -> ModularQExpansionCoefficients {
        ModularQExpansionCoefficients::j_invariant_current_table()
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
    pub fn coefficients(
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
    pub fn from_tau(
        tau: UpperHalfPlanePoint,
        truncation: QExpansionTruncation,
    ) -> Result<JInvariantQExpansionApprox, AnalyticCurveError> {
        let family = Self::new();
        <Self as ModularQExpansionFamily>::evaluate_tau(&family, tau, truncation)
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
