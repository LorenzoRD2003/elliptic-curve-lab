/// Tentative projective point representation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProjectivePoint<Scalar> {
    pub x: Scalar,
    pub y: Scalar,
    pub z: Scalar,
}

impl<Scalar> ProjectivePoint<Scalar> {
    /// Builds a projective point.
    pub fn new(x: Scalar, y: Scalar, z: Scalar) -> Self {
        Self { x, y, z }
    }
}
