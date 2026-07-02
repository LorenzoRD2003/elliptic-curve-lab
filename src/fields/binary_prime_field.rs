use crate::fields::traits::*;
use core::fmt;
use core::hash::Hash;
use core::num::NonZeroU32;

use num_bigint::{BigInt, BigUint};

use crate::fields::{
    FieldCharacteristic, FieldError,
    traits::{
        CbrtField, EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
    },
};

/// The prime field `𝔽₂`.
#[derive(Clone, Copy, Debug)]
pub struct Fp2;

/// Element of [`Fp2`], represented by `0` or `1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Fp2Elem(u8);

impl fmt::Display for Fp2Elem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (mod 2)", self.0)
    }
}

impl Fp2Elem {
    /// Returns the canonical representative in `{0, 1}`.
    pub const fn value(&self) -> u8 {
        self.0
    }
}

impl Field for Fp2 {
    const IS_ALGEBRAICALLY_CLOSED: bool = false;

    type Elem = Fp2Elem;

    fn characteristic() -> FieldCharacteristic {
        FieldCharacteristic::Positive(BigUint::from(2u8))
    }

    fn zero() -> Self::Elem {
        Fp2Elem(0)
    }

    fn one() -> Self::Elem {
        Fp2Elem(1)
    }

    fn from_bigint(n: &BigInt) -> Self::Elem {
        if (n & BigInt::from(1u8)) == BigInt::from(0u8) {
            Fp2Elem(0)
        } else {
            Fp2Elem(1)
        }
    }

    fn add(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Fp2Elem(x.0 ^ y.0)
    }

    fn sub(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Self::add(x, y)
    }

    fn mul(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Fp2Elem(x.0 & y.0)
    }

    fn neg(x: &Self::Elem) -> Self::Elem {
        *x
    }

    fn inv(x: &Self::Elem) -> Option<Self::Elem> {
        (x.0 == 1).then_some(*x)
    }

    fn eq(x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }

    fn inverse(x: &Self::Elem) -> Result<Self::Elem, FieldError> {
        Self::inv(x).ok_or(FieldError::DivisionByZero)
    }
}

impl FiniteField for Fp2 {
    fn extension_degree() -> NonZeroU32 {
        NonZeroU32::MIN
    }

    fn cardinality_biguint() -> BigUint {
        BigUint::from(2u8)
    }

    fn check_structure() -> Result<(), FieldError> {
        Ok(())
    }
}

impl EnumerableFiniteField for Fp2 {
    fn elements() -> Vec<Self::Elem> {
        vec![Self::zero(), Self::one()]
    }
}

impl SqrtField for Fp2 {
    fn sqrt(x: &Self::Elem) -> Option<Self::Elem> {
        Some(*x)
    }
}

impl CbrtField for Fp2 {
    fn cbrt(x: &Self::Elem) -> Option<Self::Elem> {
        Some(*x)
    }
}

impl QuadraticCharacterFiniteField for Fp2 {}
