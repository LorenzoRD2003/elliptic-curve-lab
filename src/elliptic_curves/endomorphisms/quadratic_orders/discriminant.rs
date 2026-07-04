use core::hash::Hash;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use crate::elliptic_curves::{
    endomorphisms::quadratic_orders::QuadraticRadicandError,
    frobenius::FrobeniusCharacteristicPolynomial,
};
use crate::numerics::{
    is_squarefree, quadratic_radicands::ImaginaryQuadraticRadicandNormalization,
};

/// Sign classification for one integral quadratic discriminant.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiscriminantSign {
    Negative,
    Zero,
    Positive,
}

/// Integral discriminant arithmetic for the first `endomorphisms` layer.
///
/// This value object models only an integer `D`, together with lightweight
/// arithmetic and classification facts that are useful before introducing
/// fuller quadratic-order or quadratic-field abstractions.
///
/// The main example is the Frobenius discriminant `D = t^2 - 4q`,
/// where `t` is the relative Frobenius trace and `q` is the finite base-field
/// cardinality.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct QuadraticDiscriminant {
    value: BigInt,
}

impl QuadraticDiscriminant {
    /// Builds a discriminant wrapper from one integer.
    ///
    /// Complexity: `Θ(1)` once the input integer has been materialized.
    pub fn new<T: Into<BigInt>>(value: T) -> Self {
        Self {
            value: value.into(),
        }
    }

    /// Builds `D = t^2 - 4q` from explicit Frobenius-side data.
    ///
    /// Complexity: `Θ(1)`
    pub fn from_frobenius_trace_and_field_order(
        trace: impl Into<BigInt>,
        field_order: impl Into<BigInt>,
    ) -> Self {
        let trace = trace.into();
        let field_order = field_order.into();
        Self::new(&trace * &trace - BigInt::from(4u8) * field_order)
    }

    /// Builds the Frobenius-side discriminant from
    /// `χ_{π_q}(T) = T^2 - tT + q`.
    ///
    /// Concretely, if `polynomial` stores `χ_{π_q}(T)`, this returns
    /// `D = t^2 - 4q`.
    ///
    /// Complexity: `Θ(1)`
    pub fn from_frobenius_characteristic_polynomial(
        polynomial: &FrobeniusCharacteristicPolynomial,
    ) -> Self {
        Self::new(polynomial.discriminant())
    }

    /// Builds the fundamental discriminant `D_K` of the imaginary quadratic field
    /// `K = ℚ(√m)` from one integer radicand `m < 0`.
    ///
    /// The current entrypoint first reduces `m = s^2 d` to its squarefree part
    /// `d < 0` and then applies the classical maximal-order rule:
    ///
    /// - if `d ≡ 1 (mod 4)`, then `D_K = d`
    /// - otherwise, `D_K = 4d`
    ///
    /// This constructs the field discriminant of the maximal order `O_K`; it
    /// does not construct a non-maximal order `O_f`.
    ///
    /// Complexity: dominated by `num-prime`.
    pub fn fundamental_from_quadratic_radicand(
        radicand: impl Into<BigInt>,
    ) -> Result<Self, QuadraticRadicandError> {
        let normalization =
            ImaginaryQuadraticRadicandNormalization::from_imaginary_radicand(radicand)
                .map_err(QuadraticRadicandError::from)?;
        Ok(Self::new(normalization.fundamental_discriminant().clone()))
    }

    /// Returns the stored integral value.
    pub fn value(&self) -> &BigInt {
        &self.value
    }

    /// Returns whether the discriminant is negative, zero, or positive.
    ///
    /// Complexity: `Θ(1)`.
    pub fn sign(&self) -> DiscriminantSign {
        if self.value.is_negative() {
            DiscriminantSign::Negative
        } else if self.value.is_zero() {
            DiscriminantSign::Zero
        } else {
            DiscriminantSign::Positive
        }
    }

    /// Returns whether `D < 0`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn is_negative(&self) -> bool {
        self.sign() == DiscriminantSign::Negative
    }

    /// Returns whether `D = 0`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn is_zero(&self) -> bool {
        self.sign() == DiscriminantSign::Zero
    }

    /// Returns whether `D > 0`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn is_positive(&self) -> bool {
        self.sign() == DiscriminantSign::Positive
    }

    /// Returns whether `D ≡ 0 (mod 4)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn is_congruent_to_0_mod_4(&self) -> bool {
        normalized_mod_4(&self.value).is_zero()
    }

    /// Returns whether `D ≡ 1 (mod 4)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn is_congruent_to_1_mod_4(&self) -> bool {
        normalized_mod_4(&self.value).is_one()
    }

    /// Returns whether `D` is a fundamental discriminant.
    ///
    /// The current implementation uses the classical quadratic-field test:
    ///
    /// - `D ≡ 1 (mod 4)` and `D` is squarefree, or
    /// - `D = 4d` with `d ≡ 2, 3 (mod 4)` and `d` squarefree.
    ///
    /// We exclude the exceptional value `D = 1`, which satisfies the raw
    /// congruence-and-squarefree pattern but does not arise from a genuine
    /// quadratic field.
    ///
    /// Complexity: dominated by `num-prime`'s implementation.
    pub fn is_fundamental(&self) -> bool {
        if self.value.is_zero() || self.value.is_one() {
            return false;
        }

        if self.is_congruent_to_1_mod_4() {
            is_squarefree(&self.value)
        } else if self.is_congruent_to_0_mod_4() {
            let quarter = &self.value / 4u8;
            let quarter_mod_four = normalized_mod_4(&quarter);
            (quarter_mod_four == BigInt::from(2u8) || quarter_mod_four == BigInt::from(3u8))
                && is_squarefree(&quarter)
        } else {
            false
        }
    }
}

fn normalized_mod_4(value: &BigInt) -> BigInt {
    ((value % 4u8) + BigInt::from(4u8)) % 4u8
}
