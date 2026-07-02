use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::rational_function_field::RationalFunction;
use crate::fields::traits::*;

/// Element of the function field `F(E)` of a short-Weierstrass curve.
pub struct ShortWeierstrassFunction<F: Field> {
    pub(super) curve: ShortWeierstrassCurve<F>,
    pub(super) a_part: RationalFunction<F>,
    pub(super) b_part: RationalFunction<F>,
}

impl<F: Field> ShortWeierstrassFunction<F> {
    pub fn curve(&self) -> &ShortWeierstrassCurve<F> {
        &self.curve
    }

    pub fn a_part(&self) -> &RationalFunction<F> {
        &self.a_part
    }

    pub fn b_part(&self) -> &RationalFunction<F> {
        &self.b_part
    }

    pub fn is_zero(&self) -> bool {
        self.a_part.is_zero() && self.b_part.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.a_part.is_one() && self.b_part.is_zero()
    }

    pub fn is_constant(&self) -> bool {
        self.b_part.is_zero() && self.a_part.is_constant()
    }
}
