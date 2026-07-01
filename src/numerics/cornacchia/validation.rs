use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::CornacchiaError;

pub(super) fn validate_coefficient_and_target(
    d: &BigUint,
    m: &BigUint,
) -> Result<(), CornacchiaError> {
    if d.is_zero() {
        return Err(CornacchiaError::ZeroCoefficient);
    }
    if m.is_zero() {
        return Err(CornacchiaError::ZeroTarget);
    }
    if m.is_one() {
        return Err(CornacchiaError::TrivialTarget);
    }
    Ok(())
}
