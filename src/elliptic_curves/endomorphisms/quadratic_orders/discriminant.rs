use core::hash::Hash;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use crate::elliptic_curves::frobenius::FrobeniusCharacteristicPolynomial;
use crate::numerics::is_squarefree;

/// Sign classification for one integral quadratic discriminant.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiscriminantSign {
    Negative,
    Zero,
    Positive,
}

/// Residue class modulo `4` for one integral quadratic discriminant.
///
/// The classical discriminant congruence conditions for quadratic orders are
/// the distinguished classes `0` and `1` modulo `4`. This enum keeps those
/// two classes explicit while still remaining honest about arbitrary integers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum QuadraticDiscriminantMod4 {
    Zero,
    One,
    Other(u8),
}

/// Integral discriminant arithmetic for the first `endomorphisms` layer.
///
/// This value object models only an integer `D`, together with lightweight
/// arithmetic and classification facts that are useful before introducing
/// fuller quadratic-order or quadratic-field abstractions.
///
/// A primary motivating example is the Frobenius-side discriminant
///
/// `D = t^2 - 4q`,
///
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
    /// Complexity: `Θ(1)` big-integer arithmetic on one trace square and one
    /// multiplication by `4`.
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
    /// Complexity: `Θ(1)` after extracting the existing big-integer
    /// polynomial discriminant.
    pub fn from_frobenius_characteristic_polynomial(
        polynomial: &FrobeniusCharacteristicPolynomial,
    ) -> Self {
        Self::new(polynomial.discriminant())
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

    /// Returns the residue class of `D` modulo `4`.
    ///
    /// Complexity: `Θ(1)` big-integer remainder arithmetic.
    pub(crate) fn mod_4_class(&self) -> QuadraticDiscriminantMod4 {
        let residue = ((&self.value % 4u8) + BigInt::from(4u8)) % 4u8;
        if residue.is_zero() {
            QuadraticDiscriminantMod4::Zero
        } else if residue.is_one() {
            QuadraticDiscriminantMod4::One
        } else {
            QuadraticDiscriminantMod4::Other(
                residue
                    .try_into()
                    .expect("a residue modulo 4 should fit into u8"),
            )
        }
    }

    /// Returns whether `D ≡ 0 (mod 4)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn is_congruent_to_0_mod_4(&self) -> bool {
        self.mod_4_class() == QuadraticDiscriminantMod4::Zero
    }

    /// Returns whether `D ≡ 1 (mod 4)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn is_congruent_to_1_mod_4(&self) -> bool {
        self.mod_4_class() == QuadraticDiscriminantMod4::One
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

        match self.mod_4_class() {
            QuadraticDiscriminantMod4::One => is_squarefree(&self.value),
            QuadraticDiscriminantMod4::Zero => {
                let quarter = &self.value / 4u8;
                let quarter_mod_four = ((&quarter % 4u8) + BigInt::from(4u8)) % 4u8;
                (quarter_mod_four == BigInt::from(2u8) || quarter_mod_four == BigInt::from(3u8))
                    && is_squarefree(&quarter)
            }
            QuadraticDiscriminantMod4::Other(_) => false,
        }
    }
}
