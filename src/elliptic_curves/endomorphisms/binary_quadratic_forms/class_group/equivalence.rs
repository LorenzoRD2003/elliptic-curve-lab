use num_bigint::BigInt;
use num_traits::One;

use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::BinaryQuadraticForm;
use crate::numerics::{extended_gcd_bigint, gcd_bigint};

/// Applies a proper unimodular change of variables to a binary quadratic form.
///
/// Given coprime integers `p` and `r`, this completes them to a matrix
/// `[[p,q],[r,s]]` with determinant `1` and returns the coefficients of
/// `f(px+qy,rx+sy)`. Proper equivalence preserves discriminant and class.
///
/// Complexity: `Θ(1)` coefficient arithmetic plus one extended gcd.
pub(in crate::elliptic_curves::endomorphisms::binary_quadratic_forms) fn properly_equivalent_form(
    form: &BinaryQuadraticForm,
    p: BigInt,
    matrix_r: BigInt,
) -> Option<BinaryQuadraticForm> {
    if gcd_bigint(&p, &matrix_r) != BigInt::one() {
        return None;
    }

    let (gcd, s, bezout_r) = extended_gcd_bigint(p.clone(), matrix_r.clone());
    debug_assert_eq!(gcd, BigInt::one());
    let q = -bezout_r;
    let (a, b, c) = form.coefficients();

    let transformed_a = a * &p * &p + b * &p * &matrix_r + c * &matrix_r * &matrix_r;
    let transformed_b = BigInt::from(2u8) * a * &p * &q
        + b * (&p * &s + &q * &matrix_r)
        + BigInt::from(2u8) * c * &matrix_r * &s;
    let transformed_c = a * &q * &q + b * &q * &s + c * &s * &s;

    Some(BinaryQuadraticForm::new(
        transformed_a,
        transformed_b,
        transformed_c,
    ))
}
