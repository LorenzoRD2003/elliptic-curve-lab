use crate::fields::errors::FieldError;

/// Returns `value mod modulus` after checking that the modulus is usable.
pub fn reduce_u128(value: u128, modulus: u64) -> Result<u64, FieldError> {
    if !is_valid_field_modulus(modulus) {
        return Err(FieldError::InvalidModulus { modulus });
    }

    Ok((value % u128::from(modulus)) as u64)
}

/// Returns whether the modulus is large enough to define a non-trivial field.
pub fn is_valid_field_modulus(modulus: u64) -> bool {
    modulus >= 2
}

/// Returns whether `value` is prime using trial division.
///
/// This helper is intentionally simple and deterministic. It is appropriate for
/// structural validation during scaffolding, not for high-performance primality
/// testing of large integers.
pub fn is_prime_u64(value: u64) -> bool {
    if value < 2 {
        return false;
    }

    if value == 2 {
        return true;
    }

    if value.is_multiple_of(2) {
        return false;
    }

    let mut divisor = 3_u64;
    while divisor <= value / divisor {
        if value.is_multiple_of(divisor) {
            return false;
        }
        divisor += 2;
    }

    true
}

/// Computes the extended greatest common divisor of `a` and `b`.
///
/// The returned triple `(g, x, y)` satisfies `a * x + b * y = g`, where `g` is
/// the non-negative greatest common divisor of `a` and `b`.
pub fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
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
