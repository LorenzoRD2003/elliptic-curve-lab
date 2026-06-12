use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::traits::GroupCurveModel;
use crate::numerics::PrimePowerTable;
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// Internal report for recovering the exact exponent in a cyclic `ℓ`-group.
///
/// This models the local situation where one already knows that
/// `[ℓ^e]Q = O`. The remaining task is to recover the exact `a` with
/// `ord(Q) = ℓ^a` by repeated multiplication by the prime `ℓ`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CyclicPrimaryOrderReport {
    prime: BigUint,
    exponent_bound: u32,
    exact_exponent: u32,
}

impl CyclicPrimaryOrderReport {
    pub(crate) fn prime(&self) -> &BigUint {
        &self.prime
    }

    pub(crate) fn exponent_bound(&self) -> u32 {
        self.exponent_bound
    }

    pub(crate) fn exact_exponent(&self) -> u32 {
        self.exact_exponent
    }

    pub(crate) fn removed_exponent(&self) -> u32 {
        self.exponent_bound - self.exact_exponent
    }
}

/// Multiplies one curve point by a non-negative `BigUint` scalar.
///
/// This is the internal big-integer analogue of the public `u64`
/// double-and-add surface on [`GroupCurveModel`]. Keeping it here lets
/// additive-group algorithms such as order recovery work over genuinely large
/// integers without forcing the public group trait to change all at once.
pub(crate) fn mul_scalar_biguint<C: GroupCurveModel + ?Sized>(
    curve: &C,
    point: &C::Point,
    scalar: &BigUint,
) -> Result<C::Point, CurveError> {
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    let mut result = curve.identity();
    let mut base = point.clone();
    let mut k = scalar.clone();

    while !k.is_zero() {
        if (&k & BigUint::one()) == BigUint::one() {
            result = curve.add(&result, &base)?;
        }

        k >>= 1usize;

        if !k.is_zero() {
            base = curve.double(&base)?;
        }
    }

    Ok(result)
}

/// Recovers the exact local exponent in a cyclic `ℓ`-group.
///
/// Precondition: `[ℓ^e]Q = O`.
///
/// The current implementation advances incrementally through the local
/// `ℓ`-power chain `Q, [ℓ]Q, [ℓ²]Q, ...`
///
/// so each step reuses the previous point instead of recomputing one larger
/// scalar multiple from the original `Q`.
pub(crate) fn recover_cyclic_primary_order<C: GroupCurveModel>(
    curve: &C,
    point: &C::Point,
    powers: &PrimePowerTable,
) -> Result<CyclicPrimaryOrderReport, CurveError> {
    if curve.is_identity(point) {
        return Ok(CyclicPrimaryOrderReport {
            prime: powers.prime().clone(),
            exponent_bound: powers.exponent_bound(),
            exact_exponent: 0,
        });
    }

    let exponent_bound = powers.exponent_bound();
    let mut current = point.clone();
    let mut exact_exponent = 0u32;

    while exact_exponent < exponent_bound {
        exact_exponent += 1;
        current = mul_scalar_biguint(curve, &current, powers.prime())?;
        if curve.is_identity(&current) {
            return Ok(CyclicPrimaryOrderReport {
                prime: powers.prime().clone(),
                exponent_bound,
                exact_exponent,
            });
        }
    }

    Err(CurveError::PointOrderMultipleDoesNotAnnihilatePoint {
        multiple: powers.power(exponent_bound).clone(),
    })
}
