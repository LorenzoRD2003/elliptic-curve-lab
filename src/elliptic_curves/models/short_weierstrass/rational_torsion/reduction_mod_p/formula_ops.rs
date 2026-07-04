use num_bigint::BigInt;

use crate::elliptic_curves::{
    CurveError,
    short_weierstrass::{
        group_law_core::ShortWeierstrassFormulaOps,
        rational_torsion::reduction_mod_p::small_prime_field::{ReductionPrime, ReductionResidue},
    },
};

pub(super) struct ReductionFormulaOps {
    prime: ReductionPrime,
}

impl ReductionFormulaOps {
    pub(super) fn new(prime: ReductionPrime) -> Self {
        Self { prime }
    }
}

impl ShortWeierstrassFormulaOps for ReductionFormulaOps {
    type Coord = ReductionResidue;

    fn add(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(self.prime.add(*left, *right))
    }

    fn sub(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(self.prime.sub(*left, *right))
    }

    fn mul(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(self.prime.mul(*left, *right))
    }

    fn inv(&self, value: &Self::Coord) -> Result<Self::Coord, CurveError> {
        self.prime
            .inv(*value)
            .ok_or(CurveError::NonInvertibleFunctionFieldElement)
    }

    fn lift_i64(&self, value: i64) -> Self::Coord {
        self.prime.reduce_bigint(&BigInt::from(value))
    }

    fn is_zero(&self, value: &Self::Coord) -> bool {
        value.is_zero()
    }

    fn eq(&self, left: &Self::Coord, right: &Self::Coord) -> bool {
        left == right
    }
}
