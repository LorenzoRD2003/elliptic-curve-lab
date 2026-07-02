use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunction;

impl<F: Field> Clone for ShortWeierstrassFunction<F> {
    fn clone(&self) -> Self {
        Self::new(
            self.curve().clone(),
            self.a_part().clone(),
            self.b_part().clone(),
        )
    }
}

impl<F: Field> PartialEq for ShortWeierstrassFunction<F> {
    fn eq(&self, other: &Self) -> bool {
        F::eq(self.curve().a(), other.curve().a())
            && F::eq(self.curve().b(), other.curve().b())
            && self.a_part() == other.a_part()
            && self.b_part() == other.b_part()
    }
}

impl<F: Field> fmt::Debug for ShortWeierstrassFunction<F>
where
    F::Elem: fmt::Debug,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ShortWeierstrassFunction")
            .field("curve_a", self.curve().a())
            .field("curve_b", self.curve().b())
            .field("b_part", &self.b_part())
            .field("a_part", &self.a_part())
            .finish()
    }
}
