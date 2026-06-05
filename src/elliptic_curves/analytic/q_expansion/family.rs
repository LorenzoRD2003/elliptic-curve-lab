use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, UpperHalfPlanePoint};

use super::{ModularQExpansionCoefficients, ModularQParameter, QExpansionTruncation};

/// Shared accessor surface for concrete modular `q`-expansion approximations.
///
/// This trait is intentionally neutral: it covers approximations coming from
/// both holomorphic modular forms such as `E₄` and `E₆`, and modular
/// functions such as `j`.
pub trait ModularQExpansionApproximation {
    /// Returns the validated upper-half-plane input and derived modular
    /// parameter used in this approximation.
    fn q_parameter(&self) -> &ModularQParameter;

    /// Returns the truncated complex value of the family at the stored `q`.
    fn value(&self) -> &Complex64;

    /// Returns the truncation policy used for this approximation.
    fn truncation(&self) -> QExpansionTruncation;

    /// Returns how many coefficients from the chosen family were retained.
    fn terms_used(&self) -> usize {
        self.truncation().terms()
    }

    /// Returns the original upper-half-plane input `τ`.
    fn tau(&self) -> &UpperHalfPlanePoint {
        self.q_parameter().tau()
    }

    /// Returns the modular parameter `q = e^{2π i τ}` used in the evaluation.
    fn q(&self) -> &Complex64 {
        self.q_parameter().q()
    }
}

/// Neutral family trait for modular `q`-expansions.
///
/// Mathematically, this is broader than “holomorphic modular form”.
///
/// For a modular transformation
/// `γ = [[a, b], [c, d]] ∈ SL₂(ℤ)` acting by
/// `γτ = (aτ + b) / (cτ + d)`,
/// a modular form of weight `k` satisfies
/// `f(γτ) = (cτ + d)^k f(τ)`.
///
/// The Eisenstein series `E₄` and `E₆` are holomorphic modular forms with
/// weights `4` and `6`:
///
/// - `E₄(γτ) = (cτ + d)^4 E₄(τ)`
/// - `E₆(γτ) = (cτ + d)^6 E₆(τ)`
///
/// By contrast, the classical modular `j`-invariant satisfies
/// `j(γτ) = j(τ)`. It is therefore a modular function of weight `0`, not a
/// holomorphic modular form. It also has a pole at the cusp `i∞`, which is
/// why its `q`-expansion begins with the principal part `q^{-1}` rather than
/// only nonnegative powers.
///
/// So this trait deliberately models a common `q`-expansion interface for
/// families such as `j(q)`, `E₄(q)`, and `E₆(q)` without pretending they are
/// all modular forms of the same kind.
pub trait ModularQExpansionFamily {
    type Approximation: ModularQExpansionApproximation;

    /// Returns the exact coefficient table used by this truncated family
    /// evaluation.
    fn coefficients(
        &self,
        truncation: QExpansionTruncation,
    ) -> Result<ModularQExpansionCoefficients, AnalyticCurveError>;

    /// Builds the family-specific approximation report from a shared
    /// evaluation result.
    fn build_approximation(
        &self,
        q_parameter: ModularQParameter,
        value: Complex64,
        truncation: QExpansionTruncation,
    ) -> Self::Approximation;

    /// Evaluates the chosen truncated family directly at the supplied
    /// modular parameter.
    fn evaluate_at_q(
        &self,
        q_parameter: &ModularQParameter,
        truncation: QExpansionTruncation,
    ) -> Result<Complex64, AnalyticCurveError> {
        self.coefficients(truncation)?.evaluate_at(*q_parameter.q())
    }

    /// Builds the family-specific approximation at a point `τ` in the upper
    /// half-plane.
    fn evaluate_tau(
        &self,
        tau: UpperHalfPlanePoint,
        truncation: QExpansionTruncation,
    ) -> Result<Self::Approximation, AnalyticCurveError> {
        let q_parameter = ModularQParameter::from_tau(tau);
        let value = self.evaluate_at_q(&q_parameter, truncation)?;

        Ok(self.build_approximation(q_parameter, value, truncation))
    }
}

macro_rules! impl_modular_q_expansion_accessors {
    ($approximation:ty) => {
        impl $approximation {
            /// Returns the modular `q`-parameter used in this approximation.
            pub fn q_parameter(&self) -> &ModularQParameter {
                &self.q_parameter
            }

            /// Returns the original upper-half-plane parameter `τ`.
            pub fn tau(&self) -> &UpperHalfPlanePoint {
                self.q_parameter.tau()
            }

            /// Returns the derived modular parameter `q = e^{2πiτ}`.
            pub fn q(&self) -> &Complex64 {
                self.q_parameter.q()
            }

            /// Returns the truncated complex value of the family.
            pub fn value(&self) -> &Complex64 {
                &self.value
            }

            /// Returns the truncation policy used in this approximation.
            pub fn truncation(&self) -> QExpansionTruncation {
                self.truncation
            }

            /// Returns how many coefficients from the chosen family were
            /// retained in the current approximation.
            pub fn terms_used(&self) -> usize {
                self.truncation.terms()
            }
        }

        impl crate::elliptic_curves::analytic::q_expansion::ModularQExpansionApproximation
            for $approximation
        {
            fn q_parameter(&self) -> &ModularQParameter {
                &self.q_parameter
            }

            fn value(&self) -> &Complex64 {
                &self.value
            }

            fn truncation(&self) -> QExpansionTruncation {
                self.truncation
            }
        }
    };
}

pub(super) use impl_modular_q_expansion_accessors;
