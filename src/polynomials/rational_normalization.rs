use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use crate::{
    fields::Q,
    numerics::{gcd_bigint, lcm_bigint},
    polynomials::{DensePolynomial, IntegerPolynomial},
};

/// Clears denominators and content from a polynomial in `Q[x]`.
///
/// The returned vector is the primitive integer coefficient vector in
/// ascending degree order. Multiplying a rational polynomial by this nonzero
/// rational scalar preserves its roots, so this is the right bridge when an
/// exact integer algorithm needs the same zero set.
///
/// Complexity: `Θ(n)` exact integer gcd/lcm operations for `n` coefficients.
pub(crate) fn primitive_integer_coefficients(polynomial: &DensePolynomial<Q>) -> Vec<BigInt> {
    let coefficients = polynomial.coefficients();
    if coefficients.is_empty() {
        return Vec::new();
    }

    let denominator_lcm = coefficients
        .iter()
        .fold(BigInt::one(), |accumulator, coefficient| {
            lcm_bigint(&accumulator, coefficient.denom())
        });

    let mut cleared = coefficients
        .iter()
        .map(|coefficient| coefficient.numer() * (&denominator_lcm / coefficient.denom()))
        .collect::<Vec<_>>();

    let content = cleared
        .iter()
        .fold(BigInt::zero(), |accumulator, coefficient| {
            gcd_bigint(&accumulator, coefficient)
        });

    if !content.is_zero() && content != BigInt::one() {
        for coefficient in &mut cleared {
            *coefficient /= &content;
        }
    }

    if let Some(leading) = cleared.last()
        && leading.is_negative()
    {
        for coefficient in &mut cleared {
            *coefficient = -coefficient.clone();
        }
    }

    cleared
}

/// Builds the primitive integer polynomial with the same zero set as a
/// polynomial in `Q[x]`.
///
/// Complexity: `Θ(n)` exact integer gcd/lcm operations for `n` coefficients.
pub(crate) fn primitive_integer_polynomial(polynomial: &DensePolynomial<Q>) -> IntegerPolynomial {
    IntegerPolynomial::new(primitive_integer_coefficients(polynomial))
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use crate::{
        fields::Q,
        polynomials::{DensePolynomial, rational_normalization::primitive_integer_coefficients},
    };

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn primitive_integer_coefficients_clear_denominators_and_content() {
        let polynomial = DensePolynomial::<Q>::new(vec![q(1, 2), q(3, 4), q(-5, 4)]);

        assert_eq!(
            primitive_integer_coefficients(&polynomial),
            vec![BigInt::from(-2), BigInt::from(-3), BigInt::from(5)]
        );
    }
}
