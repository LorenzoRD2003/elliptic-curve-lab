use crate::elliptic_curves::short_weierstrass::function_fields::{
    ShortWeierstrassFunctionField, ShortWeierstrassFunctionFieldPoint,
};
use crate::fields::traits::Field;

impl<F: Field> ShortWeierstrassFunctionField<F> {
    /// Returns the generic affine point `(x, y)` of the current curve.
    pub fn generic_point(&self) -> ShortWeierstrassFunctionFieldPoint<F> {
        ShortWeierstrassFunctionFieldPoint::Affine {
            x: self.x(),
            y: self.y(),
        }
    }
}
