use num_bigint::BigUint;

/// One non-negative representation `m = x² + d y²`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct DiagonalBinaryQuadraticRepresentation {
    x: BigUint,
    y: BigUint,
}

impl DiagonalBinaryQuadraticRepresentation {
    pub(crate) fn new(x: BigUint, y: BigUint) -> Self {
        Self { x, y }
    }

    /// Returns the non-negative `x` coordinate in `m = x² + d y²`.
    pub fn x(&self) -> &BigUint {
        &self.x
    }

    /// Returns the non-negative `y` coordinate in `m = x² + d y²`.
    pub fn y(&self) -> &BigUint {
        &self.y
    }

    /// Evaluates `x² + d y²` for this representation.
    ///
    /// Complexity: `Θ(1)` exact big-integer additions/multiplications.
    pub fn value(&self, d: &BigUint) -> BigUint {
        (&self.x * &self.x) + d * &self.y * &self.y
    }
}
