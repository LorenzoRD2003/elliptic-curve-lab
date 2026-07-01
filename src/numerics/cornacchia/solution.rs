use num_bigint::BigUint;

use crate::numerics::gcd_biguint;

/// One non-negative solution to `x² + d y² = m`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CornacchiaSolution {
    x: BigUint,
    y: BigUint,
}

impl CornacchiaSolution {
    pub(crate) fn new(x: BigUint, y: BigUint) -> Self {
        Self { x, y }
    }

    pub fn x(&self) -> &BigUint {
        &self.x
    }

    pub fn y(&self) -> &BigUint {
        &self.y
    }

    pub fn is_primitive(&self) -> bool {
        gcd_biguint(&self.x, &self.y) == BigUint::from(1u8)
    }
}
