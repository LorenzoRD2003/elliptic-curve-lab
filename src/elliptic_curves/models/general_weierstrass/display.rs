use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::models::general_weierstrass::GeneralWeierstrassCurve;

impl<F: Field> Clone for GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn clone(&self) -> Self {
        Self::from_validated_coefficients_unchecked(
            self.a1().clone(),
            self.a2().clone(),
            self.a3().clone(),
            self.a4().clone(),
            self.a6().clone(),
        )
    }
}

impl<F: Field> PartialEq for GeneralWeierstrassCurve<F>
where
    F::Elem: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.a1() == other.a1()
            && self.a2() == other.a2()
            && self.a3() == other.a3()
            && self.a4() == other.a4()
            && self.a6() == other.a6()
    }
}

impl<F: Field> Eq for GeneralWeierstrassCurve<F> where F::Elem: Eq {}

impl<F: Field> fmt::Display for GeneralWeierstrassCurve<F>
where
    F::Elem: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_equation_string())
    }
}

impl<F: Field> fmt::Debug for GeneralWeierstrassCurve<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("GeneralWeierstrassCurve")
            .field(
                "equation",
                &format_args!(
                    "y^2 + ({:?})xy + ({:?})y = x^3 + ({:?})x^2 + ({:?})x + ({:?})",
                    self.a1(),
                    self.a3(),
                    self.a2(),
                    self.a4(),
                    self.a6(),
                ),
            )
            .field("a1", self.a1())
            .field("a2", self.a2())
            .field("a3", self.a3())
            .field("a4", self.a4())
            .field("a6", self.a6())
            .finish()
    }
}
