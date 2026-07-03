use num_bigint::BigUint;
use num_rational::BigRational;
use num_traits::One;

use crate::numerics::NormalizedPrimePowerFactorization;

/// Returns the least positive integer `u` such that `uᵏq ∈ ℤ`.
///
/// If `q = n/d` is normalized and `d = Π pᵢ^eᵢ`, then the returned scale is
/// `Π pᵢ^⌈eᵢ/k⌉`. The exponent `k` must be positive; this is a developer-side
/// precondition because current callers use fixed mathematical exponents.
///
/// Complexity: dominated by factoring the denominator of `q`.
pub(crate) fn rational_denominator_power_clearance(value: &BigRational, exponent: u32) -> BigUint {
    assert!(
        exponent > 0,
        "rational denominator clearance exponent must be positive"
    );

    let denominator = value
        .denom()
        .to_biguint()
        .expect("BigRational denominators are normalized as positive integers");
    if denominator.is_one() {
        return BigUint::one();
    }

    let factorization = NormalizedPrimePowerFactorization::factor(&denominator)
        .expect("nontrivial rational denominator should admit a prime-power factorization");
    factorization
        .as_slice()
        .iter()
        .fold(BigUint::one(), |scale, (prime, denominator_exponent)| {
            scale * prime.pow(denominator_exponent.div_ceil(exponent))
        })
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::rational_denominator_power_clearance;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn denominator_power_clearance_uses_prime_power_ceilings() {
        assert_eq!(
            rational_denominator_power_clearance(&q(-1, 16), 4),
            2u8.into()
        );
        assert_eq!(
            rational_denominator_power_clearance(&q(1, 64), 6),
            2u8.into()
        );
        assert_eq!(
            rational_denominator_power_clearance(&q(1, 72), 2),
            12u8.into()
        );
    }

    #[test]
    fn denominator_power_clearance_is_one_for_integral_inputs() {
        assert_eq!(
            rational_denominator_power_clearance(&q(-7, 1), 4),
            1u8.into()
        );
    }
}
