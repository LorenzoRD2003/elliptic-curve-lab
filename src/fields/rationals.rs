use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::fields::{errors::FieldError, traits::Field};

/// The field of rational numbers `Q`.
///
/// This implementation is exact and uses arbitrary-precision integers through
/// [`BigInt`] and normalized rational values through [`BigRational`]. It is
/// intended for educational algebraic experiments, especially for examples
/// where an infinite base field is conceptually clearer than a finite field.
#[derive(Clone, Debug)]
pub struct Q;

impl Field for Q {
    type Elem = BigRational;

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

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::Q;
    use crate::fields::Field;

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
}
