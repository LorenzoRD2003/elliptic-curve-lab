use core::fmt;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::traits::Field;

impl<F: Field> Clone for ShortWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn clone(&self) -> Self {
        Self::from_validated_coefficients_unchecked(self.a().clone(), self.b().clone())
    }
}

impl<F: Field> PartialEq for ShortWeierstrassCurve<F>
where
    F::Elem: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.a() == other.a() && self.b() == other.b()
    }
}

impl<F: Field> Eq for ShortWeierstrassCurve<F> where F::Elem: Eq {}

impl<F: Field> fmt::Display for ShortWeierstrassCurve<F>
where
    F::Elem: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_equation_string())
    }
}

impl<F: Field> fmt::Debug for ShortWeierstrassCurve<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ShortWeierstrassCurve")
            .field(
                "equation",
                &format_args!("y^2 = x^3 + ({:?})x + ({:?})", self.a(), self.b()),
            )
            .field("a", self.a())
            .field("b", self.b())
            .finish()
    }
}
