use crate::elliptic_curves::{
    CurveError, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::traits::*;

impl<F: Field> ShortWeierstrassFunction<F> {
    pub(crate) fn ensure_same_curve(&self, rhs: &Self) -> Result<(), CurveError> {
        if F::eq(self.curve().a(), rhs.curve().a()) && F::eq(self.curve().b(), rhs.curve().b()) {
            Ok(())
        } else {
            Err(CurveError::IncompatibleFunctionFieldCurves)
        }
    }
}
