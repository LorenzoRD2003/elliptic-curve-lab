/// Tentative affine point representation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AffinePoint<Scalar> {
    pub x: Scalar,
    pub y: Scalar,
    pub infinity: bool,
}

impl<Scalar> AffinePoint<Scalar> {
    /// Builds a finite affine point.
    pub fn new(x: Scalar, y: Scalar) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }

    /// Builds the point at infinity.
    pub fn infinity(x: Scalar, y: Scalar) -> Self {
        Self {
            x,
            y,
            infinity: true,
        }
    }
}
