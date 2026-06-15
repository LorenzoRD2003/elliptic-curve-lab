use num_bigint::BigUint;

use crate::elliptic_curves::{CurveError, traits::BigScalarGroupCurveModel};
use crate::numerics::PrimePowerTable;

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

/// Crate-private extension for recovering the exact local exponent in a cyclic
/// `ℓ`-group.
pub(crate) trait CyclicPrimaryOrderGroupCurveModel: BigScalarGroupCurveModel {
    /// Precondition: `[ℓ^e]Q = O`.
    ///
    /// The current implementation advances incrementally through the local
    /// `ℓ`-power chain `Q, [ℓ]Q, [ℓ²]Q, ...`
    ///
    /// so each step reuses the previous point instead of recomputing one
    /// larger scalar multiple from the original `Q`.
    fn recover_cyclic_primary_order(
        &self,
        point: &Self::Point,
        powers: &PrimePowerTable,
    ) -> Result<CyclicPrimaryOrderReport, CurveError> {
        if self.is_identity(point) {
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
            current = self.mul_scalar_biguint(&current, powers.prime())?;
            if self.is_identity(&current) {
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
}

impl<T: BigScalarGroupCurveModel + ?Sized> CyclicPrimaryOrderGroupCurveModel for T {}
