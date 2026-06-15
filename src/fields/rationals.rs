use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::fields::{
    error::FieldError,
    traits::{CbrtField, Field, SqrtField},
};

/// The field of rational numbers `Q`.  
///
/// This implementation is exact and uses arbitrary-precision integers through
/// [`BigInt`] and normalized rational values through [`BigRational`]. It is
/// intended for educational algebraic experiments, especially for examples
/// where an infinite base field is conceptually clearer than a finite field.
#[derive(Clone, Debug)]
pub struct Q;

impl Field for Q {
    /// The rational numbers are not algebraically closed.
    ///
    /// For example, `x^2 - 2` has no root in `Q`.
    const IS_ALGEBRAICALLY_CLOSED: bool = false;

    type Elem = BigRational;

    fn characteristic() -> u64 {
        0
    }

    fn zero() -> Self::Elem {
        BigRational::zero()
    }

    fn one() -> Self::Elem {
        BigRational::one()
    }

    fn from_i64(n: i64) -> Self::Elem {
        BigRational::from_integer(BigInt::from(n))
    }

    fn add(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        x + y
    }

    fn sub(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        x - y
    }

    fn mul(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        x * y
    }

    fn neg(x: &Self::Elem) -> Self::Elem {
        -x
    }

    fn inv(x: &Self::Elem) -> Option<Self::Elem> {
        if x.is_zero() {
            None
        } else {
            Some(x.clone().recip())
        }
    }

    fn eq(x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }

    fn inverse(x: &Self::Elem) -> Result<Self::Elem, FieldError> {
        Self::inv(x).ok_or(FieldError::DivisionByZero)
    }

    fn elem_from_u64(value: u64) -> Self::Elem {
        BigRational::from_integer(BigInt::from(value))
    }
}

impl Q {
    fn exact_integer_square_root(value: &BigInt) -> Option<BigInt> {
        if value < &BigInt::zero() {
            return None;
        }

        if value.is_zero() {
            return Some(BigInt::zero());
        }

        let one = BigInt::one();
        let mut low = BigInt::zero();
        let mut high = BigInt::one();

        while &high * &high < *value {
            high <<= 1_u32;
        }

        while low <= high {
            let mid = (&low + &high) >> 1_u32;
            let mid_squared = &mid * &mid;

            if mid_squared == *value {
                return Some(mid);
            }

            if mid_squared < *value {
                low = mid + &one;
            } else {
                high = mid - &one;
            }
        }

        None
    }

    fn exact_integer_cube_root(value: &BigInt) -> Option<BigInt> {
        if value.is_zero() {
            return Some(BigInt::zero());
        }

        let sign = if value < &BigInt::zero() { -1 } else { 1 };
        let absolute_value = if sign < 0 { -value } else { value.clone() };
        let one = BigInt::one();
        let mut low = BigInt::zero();
        let mut high = BigInt::one();

        while &high * &high * &high < absolute_value {
            high <<= 1_u32;
        }

        while low <= high {
            let mid = (&low + &high) >> 1_u32;
            let mid_cubed = &mid * &mid * &mid;

            if mid_cubed == absolute_value {
                return Some(if sign < 0 { -mid } else { mid });
            }

            if mid_cubed < absolute_value {
                low = mid + &one;
            } else {
                high = mid - &one;
            }
        }

        None
    }
}

impl SqrtField for Q {
    /// Returns an exact rational square root when the rational is already a
    /// square in `Q`.
    ///
    /// Concretely, the current implementation succeeds only when the reduced
    /// numerator and denominator are both perfect integer squares.
    ///
    /// TODO: keep this exact semantics, but replace the current generic
    /// integer-square-root helper with a more principled bigint routine if the
    /// crate later grows broader exact-number infrastructure.
    fn sqrt(x: &Self::Elem) -> Option<Self::Elem> {
        if x < &BigRational::zero() {
            return None;
        }

        let numerator = Self::exact_integer_square_root(x.numer())?;
        let denominator = Self::exact_integer_square_root(x.denom())?;
        Some(BigRational::new(numerator, denominator))
    }
}

impl CbrtField for Q {
    /// Returns an exact rational cube root when the rational is already a
    /// cube in `Q`.
    ///
    /// Concretely, the current implementation succeeds only when the reduced
    /// numerator and denominator are both perfect integer cubes.
    ///
    /// Negative inputs are allowed: if `x = -(a/b)^3`, the returned root is
    /// `-a/b`.
    fn cbrt(x: &Self::Elem) -> Option<Self::Elem> {
        let numerator = Self::exact_integer_cube_root(x.numer())?;
        let denominator = Self::exact_integer_cube_root(x.denom())?;
        Some(BigRational::new(numerator, denominator))
    }
}

#[cfg(test)]
mod tests {
    use std::hint::black_box;

    use num_bigint::BigInt;
    use num_rational::BigRational;

    use crate::fields::{
        Q,
        traits::{CbrtField, Field, SqrtField},
    };

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn zero_one_and_integer_embeddings_are_exact() {
        assert!(Q::eq(&Q::zero(), &q(0, 1)));
        assert!(Q::eq(&Q::one(), &q(1, 1)));
        assert!(Q::eq(&Q::from_i64(-7), &q(-7, 1)));
        assert!(Q::eq(&Q::elem_from_u64(42), &q(42, 1)));
    }

    #[test]
    fn arithmetic_operations_match_exact_rational_arithmetic() {
        let x = q(2, 3);
        let y = q(5, 7);

        assert!(Q::eq(&Q::add(&x, &y), &q(29, 21)));
        assert!(Q::eq(&Q::sub(&x, &y), &q(-1, 21)));
        assert!(Q::eq(&Q::mul(&x, &y), &q(10, 21)));
        assert!(Q::eq(&Q::neg(&x), &q(-2, 3)));
    }

    #[test]
    fn inverses_and_division_behave_as_expected() {
        let x = q(3, 5);
        let y = q(7, 11);

        assert!(Q::eq(
            &Q::inv(&x).expect("non-zero rational should be invertible"),
            &q(5, 3),
        ));
        assert!(Q::eq(
            &Q::inverse(&y).expect("non-zero rational should be invertible"),
            &q(11, 7),
        ));
        assert!(Q::eq(
            &Q::div(&x, &y).expect("division by non-zero rational should work"),
            &q(33, 35),
        ));
    }

    #[test]
    fn inverse_and_division_reject_zero() {
        assert!(Q::inv(&Q::zero()).is_none());
        assert!(matches!(
            Q::inverse(&Q::zero()),
            Err(crate::fields::FieldError::DivisionByZero)
        ));
        assert!(matches!(
            Q::div(&Q::one(), &Q::zero()),
            Err(crate::fields::FieldError::DivisionByZero)
        ));
    }

    #[test]
    fn equality_uses_canonical_rational_normalization() {
        assert!(Q::eq(&q(1, 2), &q(2, 4)));
        assert!(Q::eq(&q(-3, 9), &q(-1, 3)));
        assert!(!Q::eq(&q(1, 2), &q(3, 5)));
    }

    #[test]
    fn square_cube_and_pow_are_exact() {
        let x = q(-2, 3);

        assert!(Q::eq(&Q::square(&x), &q(4, 9)));
        assert!(Q::eq(&Q::cube(&x), &q(-8, 27)));
        assert!(Q::eq(&Q::pow(&x, 0), &q(1, 1)));
        assert!(Q::eq(&Q::pow(&x, 4), &q(16, 81)));
    }

    #[test]
    fn additive_and_multiplicative_identities_hold_on_samples() {
        let samples = [q(0, 1), q(1, 2), q(-7, 5), q(13, 9)];

        for sample in samples {
            assert!(Q::eq(&Q::add(&sample, &Q::zero()), &sample));
            assert!(Q::eq(&Q::mul(&sample, &Q::one()), &sample));
        }
    }

    #[test]
    fn algebraic_closedness_metadata_matches_q() {
        assert!(!black_box(Q::IS_ALGEBRAICALLY_CLOSED));
    }

    #[test]
    fn sqrt_returns_exact_roots_for_rational_squares() {
        let root = Q::sqrt(&q(4, 9)).expect("4/9 should be a square in Q");

        assert!(Q::eq(&root, &q(2, 3)));
        assert!(Q::eq(&Q::square(&root), &q(4, 9)));
    }

    #[test]
    fn sqrt_rejects_negative_and_non_square_rationals() {
        assert!(Q::sqrt(&q(-1, 1)).is_none());
        assert!(Q::sqrt(&q(2, 1)).is_none());
        assert!(!Q::has_square_root(&q(2, 1)));
    }

    #[test]
    fn sqrt_pair_returns_opposite_rational_roots() {
        let (left, right) = Q::sqrt_pair(&q(9, 16)).expect("9/16 should be a square in Q");

        assert!(Q::eq(&Q::square(&left), &q(9, 16)));
        assert!(Q::eq(&Q::square(&right), &q(9, 16)));
        assert!(Q::eq(&right, &Q::neg(&left)));
    }

    #[test]
    fn cbrt_returns_exact_roots_for_rational_cubes() {
        let root = Q::cbrt(&q(-8, 27)).expect("-8/27 should be a cube in Q");

        assert!(Q::eq(&root, &q(-2, 3)));
        assert!(Q::eq(&Q::cube(&root), &q(-8, 27)));
    }

    #[test]
    fn cbrt_rejects_non_cube_rationals() {
        assert!(Q::cbrt(&q(2, 1)).is_none());
        assert!(Q::cbrt(&q(4, 9)).is_none());
        assert!(!Q::has_cube_root(&q(2, 1)));
    }
}
