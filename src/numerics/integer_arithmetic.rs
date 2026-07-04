use crate::numerics::{gcd_bigint, gcd_biguint};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Signed, Zero};

/// Returns `base^exponent` by exponentiation by squaring.
///
/// This small helper fills the gap between `BigInt::pow(u32)` and call sites
/// whose natural exponent is a `usize`, such as polynomial degrees.
///
/// Complexity: `Θ(log exponent)` exact integer multiplications.
pub(crate) fn pow_bigint_usize(base: &BigInt, exponent: usize) -> BigInt {
    let mut result = BigInt::one();
    let mut power = base.clone();
    let mut exponent = exponent;

    while exponent > 0 {
        if exponent % 2 == 1 {
            result *= &power;
        }
        exponent /= 2;
        if exponent > 0 {
            power = &power * &power;
        }
    }

    result
}

/// Returns the least common multiple of two nonnegative integers.
///
/// By convention this returns `0` if either input is `0`.
///
/// Complexity: one `gcd` computation plus exact integer multiplication and division.
pub(crate) fn lcm_biguint(left: &BigUint, right: &BigUint) -> BigUint {
    if left.is_zero() || right.is_zero() {
        return BigUint::zero();
    }

    (left / gcd_biguint(left, right)) * right
}

/// Returns the least common multiple of two integers.
///
/// The result is always nonnegative. By convention this returns `0` if either
/// input is `0`.
///
/// Complexity: one `gcd` computation plus exact integer multiplication and
/// division.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn lcm_bigint(left: &BigInt, right: &BigInt) -> BigInt {
    if left.is_zero() || right.is_zero() {
        return BigInt::zero();
    }

    ((left / gcd_bigint(left, right)) * right).abs()
}

/// Returns the LCM of a finite set of nonnegative integers.
///
/// The empty set returns `1`.
///
/// Complexity: `Θ(k)` calls to [`lcm_biguint`] for `k` inputs.
/// If `n = max(values)`, then it is `O(k log n)`.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn lcm_biguints<'a, I>(values: I) -> BigUint
where
    I: IntoIterator<Item = &'a BigUint>,
{
    values
        .into_iter()
        .fold(BigUint::from(1u8), |accumulator, value| {
            lcm_biguint(&accumulator, value)
        })
}

/// Returns the greatest common divisor of two nonnegative machine integers.
///
/// Complexity: `Θ(log min(a, b))` exact remainder steps.
pub(crate) fn gcd_usize(mut left: usize, mut right: usize) -> usize {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

/// Returns the least common multiple of two nonnegative machine integers.
///
/// By convention this returns `0` if either input is `0`.
///
/// Complexity: one `gcd` computation plus exact integer multiplication and
/// division.
pub(crate) fn lcm_usize(left: usize, right: usize) -> usize {
    if left == 0 || right == 0 {
        return 0;
    }

    left / gcd_usize(left, right) * right
}

/// Returns the quotients `n / p` for the distinct prime divisors `p | n`.
///
/// This is the exact family of candidates needed by the standard point-order
/// test: a point killed by `[n]` has exact order `n` exactly when it is not
/// already killed by any quotient `n / p` with prime `p | n`.
///
/// Examples:
/// - `n = 2` returns `[1]`
/// - `n = 12 = 2^2 * 3` returns `[6, 4]`
/// - `n = 27 = 3^3` returns `[9]`
///
/// Complexity: `Θ(sqrt(n))` divisibility checks in the current trial-division
/// implementation.
pub(crate) fn quotients_by_distinct_prime_factors(n: usize) -> Vec<usize> {
    let mut quotients = Vec::new();
    let mut remaining = n;
    let mut factor = 2;

    while factor * factor <= remaining {
        if remaining.is_multiple_of(factor) {
            quotients.push(n / factor);
            while remaining.is_multiple_of(factor) {
                remaining /= factor;
            }
        }

        factor += if factor == 2 { 1 } else { 2 };
    }

    if remaining > 1 && n > 1 {
        quotients.push(n / remaining);
    }

    quotients
}

#[cfg(test)]
mod tests {

    use super::{
        gcd_usize, lcm_bigint, lcm_biguint, lcm_biguints, lcm_usize, pow_bigint_usize,
        quotients_by_distinct_prime_factors,
    };
    use crate::numerics::gcd_biguint;
    use num_bigint::{BigInt, BigUint};
    use num_traits::One;

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    #[test]
    fn gcd_biguint_handles_zero_and_nontrivial_inputs() {
        assert_eq!(gcd_biguint(&bu(0), &bu(0)), bu(0));
        assert_eq!(gcd_biguint(&bu(0), &bu(42)), bu(42));
        assert_eq!(gcd_biguint(&bu(42), &bu(0)), bu(42));
        assert_eq!(gcd_biguint(&bu(84), &bu(126)), bu(42));
    }

    #[test]
    fn lcm_biguint_uses_the_zero_convention_and_exact_arithmetic() {
        assert_eq!(lcm_biguint(&bu(0), &bu(42)), bu(0));
        assert_eq!(lcm_biguint(&bu(12), &bu(18)), bu(36));
        assert_eq!(lcm_biguint(&bu(8), &bu(9)), bu(72));
    }

    #[test]
    fn lcm_bigint_normalizes_signs() {
        assert_eq!(
            lcm_bigint(&BigInt::from(-21), &BigInt::from(6)),
            BigInt::from(42)
        );
        assert_eq!(
            lcm_bigint(&BigInt::from(0), &BigInt::from(-6)),
            BigInt::from(0)
        );
    }

    #[test]
    fn lcm_biguints_accumulates_over_a_family() {
        let values = [bu(6), bu(10), bu(15)];

        assert_eq!(lcm_biguints(values.iter()), bu(30));
    }

    #[test]
    fn lcm_biguints_returns_one_on_the_empty_family() {
        let values = Vec::<BigUint>::new();

        assert_eq!(lcm_biguints(values.iter()), bu(1));
    }

    #[test]
    fn gcd_and_lcm_usize_follow_the_expected_small_integer_conventions() {
        assert_eq!(gcd_usize(84, 126), 42);
        assert_eq!(lcm_usize(0, 42), 0);
        assert_eq!(lcm_usize(12, 18), 36);
    }

    #[test]
    fn quotients_by_distinct_prime_factors_match_the_exact_order_candidates() {
        assert_eq!(quotients_by_distinct_prime_factors(2), vec![1]);
        assert_eq!(quotients_by_distinct_prime_factors(12), vec![6, 4]);
        assert_eq!(quotients_by_distinct_prime_factors(27), vec![9]);
    }

    #[test]
    fn pow_bigint_usize_uses_exponentiation_by_squaring() {
        assert_eq!(pow_bigint_usize(&BigInt::from(-2), 0), BigInt::one());
        assert_eq!(pow_bigint_usize(&BigInt::from(-2), 5), BigInt::from(-32));
        assert_eq!(pow_bigint_usize(&BigInt::from(3), 10), BigInt::from(59_049));
    }
}
