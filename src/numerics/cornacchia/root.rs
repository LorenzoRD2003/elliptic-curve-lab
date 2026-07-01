use num_bigint::BigUint;
use num_traits::Zero;

use super::CornacchiaError;

pub(super) fn validate_root(d: &BigUint, m: &BigUint, r: &BigUint) -> Result<(), CornacchiaError> {
    if ((r % m) * (r % m) + (d % m)) % m != BigUint::zero() {
        return Err(CornacchiaError::RootDoesNotSolveCongruence);
    }
    Ok(())
}

pub(super) fn normalize_root(r: &BigUint, m: &BigUint) -> BigUint {
    let normalized = r % m;
    if (&normalized * BigUint::from(2u8)) > *m {
        m - normalized
    } else {
        normalized
    }
}
