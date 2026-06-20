use crate::numerics::gcd_biguint;
use num_bigint::BigUint;
use num_traits::Zero;

/// Computes the extended greatest common divisor of `a` and `b`.
///
/// The returned triple `(g, x, y)` satisfies `a * x + b * y = g`, where `g` is
/// the non-negative greatest common divisor of `a` and `b`.
pub(crate) fn extended_gcd_i128(a: i128, b: i128) -> (i128, i128, i128) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (1_i128, 0_i128);
    let (mut old_t, mut t) = (0_i128, 1_i128);

    while r != 0 {
        let quotient = old_r / r;

        let next_r = old_r - quotient * r;
        old_r = r;
        r = next_r;

        let next_s = old_s - quotient * s;
        old_s = s;
        s = next_s;

        let next_t = old_t - quotient * t;
        old_t = t;
        t = next_t;
    }

    (old_r.abs(), old_s, old_t)
}

/// Returns the least common multiple of two nonnegative integers.
///
/// By convention this returns `0` if either input is `0`.
///
/// Complexity: one `gcd` computation plus exact integer multiplication and division.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn lcm_biguint(left: &BigUint, right: &BigUint) -> BigUint {
    if left.is_zero() || right.is_zero() {
        return BigUint::zero();
    }

    (left / gcd_biguint(left, right)) * right
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

/// Returns `value^2` as `usize`, failing if the educational exact result does
/// not fit.
pub(crate) fn square_u64_as_usize(value: u64) -> usize {
    usize::try_from(u128::from(value) * u128::from(value))
        .expect("educational exact square should fit into usize")
}

/// Returns `base^exponent` as `usize`, failing if the educational exact result
/// does not fit.
pub(crate) fn pow_u64_as_usize(base: u64, exponent: u32) -> usize {
    usize::try_from(u128::from(base).pow(exponent))
        .expect("educational exact power should fit into usize")
}

#[cfg(test)]
mod tests {
    use super::{
        gcd_usize, lcm_biguint, lcm_biguints, lcm_usize, pow_u64_as_usize,
        quotients_by_distinct_prime_factors, square_u64_as_usize,
    };
    use crate::numerics::gcd_biguint;
    use num_bigint::BigUint;

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
    fn checked_small_integer_power_helpers_preserve_exact_values() {
        assert_eq!(square_u64_as_usize(7), 49);
        assert_eq!(pow_u64_as_usize(5, 4), 625);
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
}
