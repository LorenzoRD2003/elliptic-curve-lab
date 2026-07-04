use core::fmt;

use num_bigint::{BigInt, BigUint, Sign};
use num_prime::nt_funcs::factorize;
use num_traits::{One, Zero};

/// Integral-basis shape for the maximal order of `ℚ(√d)` with squarefree `d`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum QuadraticIntegerBasisKind {
    /// `O_K = ℤ[√d]`, equivalently the discriminant is `4d`.
    ZSqrtD,
    /// `O_K = ℤ[(1 + √d) / 2]`, equivalently the discriminant is `d`.
    ZHalfOnePlusSqrtD,
}

impl QuadraticIntegerBasisKind {
    /// Returns the maximal-order basis shape determined by one squarefree part `d`.
    ///
    /// For squarefree `d`, the classical rule is:
    /// - `O_K = ℤ[(1 + √d)/2]` when `d ≡ 1 (mod 4)`
    /// - `O_K = ℤ[√d]` otherwise
    pub(crate) fn for_squarefree_part(squarefree_part: &BigInt) -> Self {
        let normalized_mod_4 = ((squarefree_part % 4u8) + BigInt::from(4u8)) % 4u8;
        if normalized_mod_4 == BigInt::one() {
            Self::ZHalfOnePlusSqrtD
        } else {
            Self::ZSqrtD
        }
    }

    /// Returns the corresponding fundamental discriminant `D_K`.
    pub(crate) fn fundamental_discriminant_for(self, squarefree_part: &BigInt) -> BigInt {
        match self {
            Self::ZSqrtD => BigInt::from(4u8) * squarefree_part,
            Self::ZHalfOnePlusSqrtD => squarefree_part.clone(),
        }
    }

    /// Adjusts the square factor in `4s²d = f²D_K`.
    pub(crate) fn discriminant_square_root_factor(self, square_root_factor: BigUint) -> BigUint {
        match self {
            Self::ZSqrtD => square_root_factor,
            Self::ZHalfOnePlusSqrtD => square_root_factor * BigUint::from(2u8),
        }
    }
}

/// Failure modes for the current imaginary-radicand normalization helper.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ImaginaryQuadraticRadicandError {
    ZeroRadicand,
    NonNegativeRadicand,
}

impl fmt::Display for ImaginaryQuadraticRadicandError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroRadicand => {
                write!(f, "quadratic-radicand normalization is undefined for m = 0")
            }
            Self::NonNegativeRadicand => write!(
                f,
                "the current quadratic-radicand normalization helper only supports imaginary inputs m < 0"
            ),
        }
    }
}

impl std::error::Error for ImaginaryQuadraticRadicandError {}

/// Canonical normalization data for one imaginary quadratic radicand `m < 0`.
///
/// This helper records `m = s²d` with `s ≥ 1` and `d < 0` squarefree,
/// together with the `mod 4` branch that determines the maximal-order
/// discriminant of `ℚ(√m) = ℚ(√d)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ImaginaryQuadraticRadicandNormalization {
    original_radicand: BigInt,
    square_root_factor: BigUint,
    squarefree_part: BigInt,
    fundamental_discriminant: BigInt,
    integral_basis_kind: QuadraticIntegerBasisKind,
}

impl ImaginaryQuadraticRadicandNormalization {
    /// Normalizes one imaginary quadratic radicand `m < 0`.
    ///
    /// The output stores:
    /// - the decomposition `m = s²d` with squarefree `d < 0`
    /// - the classical discriminant branch
    /// - the corresponding maximal-order basis shape
    ///
    /// Complexity: dominated by `num-prime`.
    pub(crate) fn from_imaginary_radicand(
        radicand: impl Into<BigInt>,
    ) -> Result<Self, ImaginaryQuadraticRadicandError> {
        let original_radicand = radicand.into();
        Self::validate_imaginary_radicand(&original_radicand)?;

        let (square_root_factor, squarefree_magnitude) =
            split_square_part(original_radicand.magnitude());
        let squarefree_part = -BigInt::from(squarefree_magnitude);

        let integral_basis_kind = QuadraticIntegerBasisKind::for_squarefree_part(&squarefree_part);
        let fundamental_discriminant =
            integral_basis_kind.fundamental_discriminant_for(&squarefree_part);

        Ok(Self {
            original_radicand,
            square_root_factor,
            squarefree_part,
            fundamental_discriminant,
            integral_basis_kind,
        })
    }

    /// Returns the fundamental discriminant `D_K` of `Q(√m) = Q(√d)`.
    pub(crate) fn fundamental_discriminant(&self) -> &BigInt {
        &self.fundamental_discriminant
    }

    fn validate_imaginary_radicand(
        radicand: &BigInt,
    ) -> Result<(), ImaginaryQuadraticRadicandError> {
        match radicand.sign() {
            Sign::NoSign => Err(ImaginaryQuadraticRadicandError::ZeroRadicand),
            Sign::Plus => Err(ImaginaryQuadraticRadicandError::NonNegativeRadicand),
            Sign::Minus => Ok(()),
        }
    }
}

/// Splits one non-negative integer as `n = s²d` with squarefree `d`.
///
/// The returned pair is `(s, d)`, where `s ≥ 1` and `d ≥ 0` is squarefree.
/// By convention, `0` is represented as `0 = 1² * 0`.
///
/// Complexity: dominated by `num-prime`.
pub(crate) fn split_square_part(value: &BigUint) -> (BigUint, BigUint) {
    if value.is_zero() {
        return (BigUint::one(), BigUint::zero());
    }

    let factors = factorize(value.clone());
    let mut square_root_factor = BigUint::one();
    let mut squarefree_part = BigUint::one();

    for (prime, exponent) in factors {
        let half = exponent / 2;
        if half > 0 {
            square_root_factor *= prime.pow(half as u32);
        }
        if exponent % 2 == 1 {
            squarefree_part *= prime;
        }
    }
    (square_root_factor, squarefree_part)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_square_part_keeps_zero_total_by_convention() {
        assert_eq!(
            split_square_part(&BigUint::from(0u8)),
            (BigUint::from(1u8), BigUint::from(0u8))
        );
    }

    #[test]
    fn split_square_part_extracts_the_maximal_square_factor() {
        assert_eq!(
            split_square_part(&BigUint::from(108u16)),
            (BigUint::from(6u8), BigUint::from(3u8))
        );
        assert_eq!(
            split_square_part(&BigUint::from(45u8)),
            (BigUint::from(3u8), BigUint::from(5u8))
        );
    }

    #[test]
    fn basis_kind_selection_follows_the_classical_mod_four_rule() {
        assert_eq!(
            QuadraticIntegerBasisKind::for_squarefree_part(&BigInt::from(-3)),
            QuadraticIntegerBasisKind::ZHalfOnePlusSqrtD
        );
        assert_eq!(
            QuadraticIntegerBasisKind::for_squarefree_part(&BigInt::from(-5)),
            QuadraticIntegerBasisKind::ZSqrtD
        );
    }

    #[test]
    fn basis_kind_recovers_the_matching_fundamental_discriminant() {
        assert_eq!(
            QuadraticIntegerBasisKind::ZHalfOnePlusSqrtD
                .fundamental_discriminant_for(&BigInt::from(-3)),
            BigInt::from(-3)
        );
        assert_eq!(
            QuadraticIntegerBasisKind::ZSqrtD.fundamental_discriminant_for(&BigInt::from(-5)),
            BigInt::from(-20)
        );
    }

    #[test]
    fn basis_kind_adjusts_the_discriminant_square_factor() {
        assert_eq!(
            QuadraticIntegerBasisKind::ZSqrtD.discriminant_square_root_factor(BigUint::from(3u8)),
            BigUint::from(3u8)
        );
        assert_eq!(
            QuadraticIntegerBasisKind::ZHalfOnePlusSqrtD
                .discriminant_square_root_factor(BigUint::from(3u8)),
            BigUint::from(6u8)
        );
    }

    #[test]
    fn normalization_rejects_zero_and_positive_radicands() {
        assert_eq!(
            ImaginaryQuadraticRadicandNormalization::from_imaginary_radicand(0),
            Err(ImaginaryQuadraticRadicandError::ZeroRadicand)
        );
        assert_eq!(
            ImaginaryQuadraticRadicandNormalization::from_imaginary_radicand(5),
            Err(ImaginaryQuadraticRadicandError::NonNegativeRadicand)
        );
    }

    #[test]
    fn normalization_of_minus_one_uses_the_z_sqrt_d_branch() {
        let normalization = ImaginaryQuadraticRadicandNormalization::from_imaginary_radicand(-1)
            .expect("-1 should define an imaginary quadratic field");

        assert_eq!(normalization.original_radicand, BigInt::from(-1));
        assert_eq!(normalization.square_root_factor, BigUint::from(1u8));
        assert_eq!(normalization.squarefree_part, BigInt::from(-1));
        assert_eq!(normalization.fundamental_discriminant(), &BigInt::from(-4));
        assert_eq!(
            normalization.integral_basis_kind,
            QuadraticIntegerBasisKind::ZSqrtD
        );
    }

    #[test]
    fn normalization_of_minus_three_uses_the_half_integral_branch() {
        let normalization = ImaginaryQuadraticRadicandNormalization::from_imaginary_radicand(-3)
            .expect("-3 should define an imaginary quadratic field");

        assert_eq!(normalization.square_root_factor, BigUint::from(1u8));
        assert_eq!(normalization.squarefree_part, BigInt::from(-3));
        assert_eq!(normalization.fundamental_discriminant(), &BigInt::from(-3));
        assert_eq!(
            normalization.integral_basis_kind,
            QuadraticIntegerBasisKind::ZHalfOnePlusSqrtD
        );
    }

    #[test]
    fn normalization_extracts_the_square_part_before_applying_the_mod_four_rule() {
        let normalization = ImaginaryQuadraticRadicandNormalization::from_imaginary_radicand(-12)
            .expect("-12 should normalize through the squarefree part -3");

        assert_eq!(normalization.square_root_factor, BigUint::from(2u8));
        assert_eq!(normalization.squarefree_part, BigInt::from(-3));
        assert_eq!(normalization.fundamental_discriminant(), &BigInt::from(-3));
        assert_eq!(
            normalization.integral_basis_kind,
            QuadraticIntegerBasisKind::ZHalfOnePlusSqrtD
        );
    }

    #[test]
    fn normalization_handles_the_four_d_branch_after_squarefree_reduction() {
        let normalization = ImaginaryQuadraticRadicandNormalization::from_imaginary_radicand(-20)
            .expect("-20 should normalize through the squarefree part -5");

        assert_eq!(normalization.square_root_factor, BigUint::from(2u8));
        assert_eq!(normalization.squarefree_part, BigInt::from(-5));
        assert_eq!(normalization.fundamental_discriminant(), &BigInt::from(-20));
        assert_eq!(
            normalization.integral_basis_kind,
            QuadraticIntegerBasisKind::ZSqrtD
        );
    }

    #[test]
    fn normalization_reconstructs_the_original_radicand() {
        let normalization = ImaginaryQuadraticRadicandNormalization::from_imaginary_radicand(-108)
            .expect("-108 should define an imaginary quadratic field");
        let square = BigInt::from(normalization.square_root_factor.pow(2));

        assert_eq!(
            square * normalization.squarefree_part,
            normalization.original_radicand
        );
    }
}
