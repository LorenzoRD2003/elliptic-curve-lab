use num_bigint::BigUint;

use crate::elliptic_curves::ShortWeierstrassCurve;

pub(super) type F7 = crate::fields::Fp7;

pub(super) fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

pub(super) fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
}
