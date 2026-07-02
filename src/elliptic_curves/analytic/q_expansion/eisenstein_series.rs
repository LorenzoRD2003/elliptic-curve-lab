use num_bigint::BigInt;
use num_complex::Complex64;
use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, UpperHalfPlanePoint,
    q_expansion::{
        ModularQExpansionCoefficients, ModularQExpansionFamily, ModularQParameter,
        QExpansionTruncation, family::impl_modular_q_expansion_accessors,
    },
};
use crate::numerics::{bernoulli_number, sigma_power_sums_up_to};

/// Validated classical holomorphic Eisenstein weight `k`.
///
/// - `k` must be even
/// - `k >= 4`
///
/// So the quasimodular exceptional case `E₂` is intentionally excluded from
/// this first holomorphic `q`-expansion family.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EisensteinSeriesWeight {
    k: u32,
}

impl EisensteinSeriesWeight {
    /// Builds a validated classical holomorphic Eisenstein weight.
    pub fn new(k: u32) -> Result<Self, AnalyticCurveError> {
        if k >= 4 && k.is_multiple_of(2) {
            Ok(Self { k })
        } else {
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        }
    }

    /// Returns the classical weight `k`.
    pub fn value(self) -> u32 {
        self.k
    }

    /// Returns the standard weight-`4` choice.
    pub fn e4() -> Self {
        Self { k: 4 }
    }

    /// Returns the standard weight-`6` choice.
    pub fn e6() -> Self {
        Self { k: 6 }
    }
}

/// Runtime family object for the holomorphic Eisenstein series `E_k(q)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EisensteinSeriesQExpansion {
    weight: EisensteinSeriesWeight,
}

impl EisensteinSeriesQExpansion {
    /// Builds the `E_k(q)` family from a validated weight.
    pub fn new(weight: EisensteinSeriesWeight) -> Self {
        Self { weight }
    }

    /// Returns the standard `E₄(q)` family.
    pub fn e4() -> Self {
        Self::new(EisensteinSeriesWeight::e4())
    }

    /// Returns the standard `E₆(q)` family.
    pub fn e6() -> Self {
        Self::new(EisensteinSeriesWeight::e6())
    }

    /// Returns the validated weight carried by this family.
    pub fn weight(&self) -> EisensteinSeriesWeight {
        self.weight
    }

    /// Returns the truncated coefficient table used by this `E_k(q)` family.
    #[allow(dead_code)]
    pub(crate) fn coefficients(
        &self,
        truncation: QExpansionTruncation,
    ) -> Result<ModularQExpansionCoefficients, AnalyticCurveError> {
        <Self as ModularQExpansionFamily>::coefficients(self, truncation)
    }

    /// Evaluates the truncated `E_k(q)` approximation at an upper-half-plane point `τ`.
    pub fn evaluate_tau(
        &self,
        tau: UpperHalfPlanePoint,
        truncation: QExpansionTruncation,
    ) -> Result<EisensteinSeriesQExpansionApprox, AnalyticCurveError> {
        <Self as ModularQExpansionFamily>::evaluate_tau(self, tau, truncation)
    }
}

/// Approximation report for the truncated `q`-expansion of a holomorphic
/// Eisenstein series `E_k`.
#[derive(Clone, Debug, PartialEq)]
pub struct EisensteinSeriesQExpansionApprox {
    weight: EisensteinSeriesWeight,
    q_parameter: ModularQParameter,
    value: Complex64,
    truncation: QExpansionTruncation,
}

impl EisensteinSeriesQExpansionApprox {
    /// Returns the validated Eisenstein weight.
    pub fn weight(&self) -> EisensteinSeriesWeight {
        self.weight
    }

    /// Returns the classical integer weight `k`.
    pub fn k(&self) -> u32 {
        self.weight.value()
    }
}

impl ModularQExpansionFamily for EisensteinSeriesQExpansion {
    type Approximation = EisensteinSeriesQExpansionApprox;

    fn coefficients(
        &self,
        truncation: QExpansionTruncation,
    ) -> Result<ModularQExpansionCoefficients, AnalyticCurveError> {
        truncated_eisenstein_series_coefficients_from_weight(self.weight.value(), truncation)
    }

    fn build_approximation(
        &self,
        q_parameter: ModularQParameter,
        value: Complex64,
        truncation: QExpansionTruncation,
    ) -> Self::Approximation {
        EisensteinSeriesQExpansionApprox {
            weight: self.weight,
            q_parameter,
            value,
            truncation,
        }
    }
}

impl_modular_q_expansion_accessors!(EisensteinSeriesQExpansionApprox);

/// Builds the first truncated coefficient table of a holomorphic Eisenstein
/// series written in the form
/// `1 + scale * Σ_{n≥1} σ_power_sum(n, sigma_power) q^n`.
///
/// Here `σ_r(n) = Σ_{d|n} d^r` is the classical divisor-power sum. The table
/// starts at exponent `0` because `E_k` with even `k ≥ 4` are holomorphic at
/// the cusp, so their `q`-expansions begin with the constant term `1` rather
/// than a principal part such as `q^{-1}`.
///
/// Because this routine needs every value
/// `σ_power_sum(1, sigma_power), ..., σ_power_sum(N, sigma_power)`
/// up to the truncation bound, it intentionally uses the batched divisor-sieve
/// helper from `numerics` instead of recomputing each sigma value separately.
fn truncated_eisenstein_series_coefficients(
    sigma_power: u32,
    scale: BigRational,
    truncation: QExpansionTruncation,
) -> ModularQExpansionCoefficients {
    let mut coefficients = Vec::with_capacity(truncation.terms());
    coefficients.push(BigRational::one());
    let sigma_values = sigma_power_sums_up_to(truncation.terms().saturating_sub(1), sigma_power);

    for sigma_value in sigma_values.iter().take(truncation.terms()).skip(1) {
        let sigma = BigRational::from_integer(sigma_value.clone());
        coefficients.push(scale.clone() * sigma);
    }

    ModularQExpansionCoefficients::new(0, coefficients)
}

/// Returns the exact rational scaling factor `-2k / B_k` appearing in the
/// classical normalized Eisenstein-series expansion (for even weights `k ≥ 4`)
///
/// `E_k(q) = 1 - (2k / B_k) Σ_{n≥1} σ_{k-1}(n) q^n`
///
/// - rejects odd weights
/// - rejects `k < 4`
/// - rejects any exceptional case where `B_k = 0`
///
/// Complexity: the dominant work is the Bernoulli-number computation,
/// so this helper is `Θ(k²)` exact rational arithmetic updates.
fn eisenstein_scale(weight: u32) -> Result<BigRational, AnalyticCurveError> {
    if weight < 4 || weight % 2 == 1 {
        return Err(AnalyticCurveError::InvalidEisensteinWeight);
    }

    let b_k = bernoulli_number(weight as usize);
    if b_k.is_zero() {
        return Err(AnalyticCurveError::InvalidEisensteinWeight);
    }

    let signed_weight = i64::from(weight);
    Ok(BigRational::from_integer(BigInt::from(-2 * signed_weight)) / b_k)
}

/// Builds the first truncated coefficient table of the holomorphic Eisenstein
/// series `E_k(q)` for an even weight `k ≥ 4`.
///
/// Internally this derives the exact scaling factor from the Bernoulli number
/// `B_k`, then performs the divisor-power-sum accumulation with
/// `σ_{k-1}(n) = Σ_{d|n} d^{k-1}`.
fn truncated_eisenstein_series_coefficients_from_weight(
    weight: u32,
    truncation: QExpansionTruncation,
) -> Result<ModularQExpansionCoefficients, AnalyticCurveError> {
    Ok(truncated_eisenstein_series_coefficients(
        weight - 1,
        eisenstein_scale(weight)?,
        truncation,
    ))
}
#[cfg(test)]
mod tests {

    use num_bigint::BigInt;
    use num_rational::BigRational;

    use crate::elliptic_curves::analytic::q_expansion::QExpansionTruncation;
    use crate::elliptic_curves::analytic::q_expansion::eisenstein_series::{
        eisenstein_scale, truncated_eisenstein_series_coefficients_from_weight,
    };
    use crate::elliptic_curves::analytic::q_expansion::{
        EisensteinSeriesQExpansion, EisensteinSeriesWeight,
    };
    use crate::elliptic_curves::analytic::{AnalyticCurveError, UpperHalfPlanePoint};

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn eisenstein_scale_matches_the_classical_e4_and_e6_constants() {
        assert_eq!(eisenstein_scale(4).unwrap(), q(240, 1));
        assert_eq!(eisenstein_scale(6).unwrap(), q(-504, 1));
    }

    #[test]
    fn eisenstein_scale_rejects_invalid_weights() {
        assert_eq!(
            eisenstein_scale(0),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
        assert_eq!(
            eisenstein_scale(2),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
        assert_eq!(
            eisenstein_scale(3),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
        assert_eq!(
            eisenstein_scale(5),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
    }

    #[test]
    fn weight_driven_coefficient_table_matches_the_known_small_examples() {
        let e4 = truncated_eisenstein_series_coefficients_from_weight(
            4,
            QExpansionTruncation::new(5).unwrap(),
        )
        .unwrap();
        let e6 = truncated_eisenstein_series_coefficients_from_weight(
            6,
            QExpansionTruncation::new(5).unwrap(),
        )
        .unwrap();

        assert_eq!(
            e4.coefficients(),
            &[q(1, 1), q(240, 1), q(2160, 1), q(6720, 1), q(17_520, 1)]
        );
        assert_eq!(
            e6.coefficients(),
            &[
                q(1, 1),
                q(-504, 1),
                q(-16_632, 1),
                q(-122_976, 1),
                q(-532_728, 1),
            ]
        );
    }

    #[test]
    fn weight_twelve_keeps_the_expected_rational_coefficient_factor() {
        let coefficients = truncated_eisenstein_series_coefficients_from_weight(
            12,
            QExpansionTruncation::new(2).unwrap(),
        )
        .unwrap();

        assert_eq!(coefficients.coefficient_of(1).unwrap(), q(65_520, 691));
    }

    #[test]
    fn weight_value_object_validates_even_holomorphic_weights() {
        assert_eq!(
            EisensteinSeriesWeight::new(2),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
        assert_eq!(
            EisensteinSeriesWeight::new(5),
            Err(AnalyticCurveError::InvalidEisensteinWeight)
        );
        assert_eq!(EisensteinSeriesWeight::new(4).unwrap().value(), 4);
        assert_eq!(EisensteinSeriesWeight::new(6).unwrap().value(), 6);
    }

    #[test]
    fn eisenstein_family_helpers_build_expected_weights() {
        assert_eq!(EisensteinSeriesQExpansion::e4().weight().value(), 4);
        assert_eq!(EisensteinSeriesQExpansion::e6().weight().value(), 6);
    }

    #[test]
    fn eisenstein_approximation_keeps_the_weight_metadata() {
        let approximation = EisensteinSeriesQExpansion::e6()
            .evaluate_tau(
                UpperHalfPlanePoint::tau_i(),
                QExpansionTruncation::new(4).unwrap(),
            )
            .unwrap();

        assert_eq!(approximation.weight(), EisensteinSeriesWeight::e6());
        assert_eq!(approximation.k(), 6);
        assert_eq!(approximation.terms_used(), 4);
    }
}
