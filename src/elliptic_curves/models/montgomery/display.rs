use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::models::montgomery::MontgomeryCurve;

impl<F: Field> Clone for MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    fn clone(&self) -> Self {
        Self::from_validated_coefficients_unchecked(self.a().clone(), self.b().clone())
    }
}

impl<F: Field> PartialEq for MontgomeryCurve<F>
where
    F::Elem: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.a() == other.a() && self.b() == other.b()
    }
}

impl<F: Field> Eq for MontgomeryCurve<F> where F::Elem: Eq {}

impl<F: Field> fmt::Display for MontgomeryCurve<F>
where
    F::Elem: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_equation_string())
    }
}

impl<F: Field> fmt::Debug for MontgomeryCurve<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MontgomeryCurve")
            .field(
                "equation",
                &format_args!("({:?})y^2 = x^3 + ({:?})x^2 + x", self.b(), self.a()),
            )
            .field("a", self.a())
            .field("b", self.b())
            .finish()
    }
}
