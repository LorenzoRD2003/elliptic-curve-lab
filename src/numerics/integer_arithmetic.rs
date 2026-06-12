use num_bigint::BigUint;
use num_traits::Zero;

/// Returns the greatest common divisor of two nonnegative integers.
///
/// Complexity: `Θ(log min(a, b))` exact remainder steps.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn gcd_biguint(left: &BigUint, right: &BigUint) -> BigUint {
    let mut left = left.clone();
    let mut right = right.clone();

    while !right.is_zero() {
        let remainder = &left % &right;
        left = right;
        right = remainder;
    }

    left
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

#[cfg(test)]
mod tests {
    use super::{gcd_biguint, lcm_biguint, lcm_biguints};
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
}
