use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::{
    rational_function_field::RationalFunction,
    traits::{FiniteField, PthRootExtraction},
};

impl<F: FiniteField> PthRootExtraction for ShortWeierstrassFunction<F>
where
    RationalFunction<F>: PthRootExtraction,
{
    fn pth_root(&self) -> Option<Self> {
        let rhs_factor = self.frobenius_quadratic_rhs_factor();

        let a_root = self.a_part().pth_root()?;
        let adjusted_b = self.b_part().div(&rhs_factor).ok()?;
        let b_root = adjusted_b.pth_root()?;

        Some(Self::new(self.curve().clone(), a_root, b_root))
    }
}

impl<F: FiniteField> ShortWeierstrassFunction<F> {
    /// Returns the rational-function factor by which the `y`-part twists under
    /// inverse Frobenius on `F(E) = F(x) ⊕ yF(x)`.
    ///
    /// If `y^2 = f(x)`, then after taking a `p`-th root of
    /// `A(x) + y B(x)` one must divide the `B`-part by
    /// $f(x)^{(p-1)/2}$, because
    /// $y^p = y \cdot (y^2)^{(p-1)/2} = y \cdot f(x)^{(p-1)/2}$.
    fn frobenius_quadratic_rhs_factor(&self) -> RationalFunction<F> {
        let p = F::characteristic();
        self.curve()
            .function_field_curve_rhs_rational_function()
            .pow_u64((p - 1) / 2)
    }
}

impl<F: FiniteField> ShortWeierstrassFunction<F>
where
    F::Elem: PartialEq,
{
    /// Tries to invert one absolute-Frobenius pullback along the first
    /// Frobenius twist of the underlying short-Weierstrass curve.
    ///
    /// The `a`-part is pulled back coefficientwise, while the `b`-part must
    /// first be divided by the factor described in
    /// [`Self::frobenius_quadratic_rhs_factor`].
    pub(crate) fn inverse_absolute_frobenius_pullback_to_twist(
        &self,
        codomain_twist: &ShortWeierstrassCurve<F>,
    ) -> Option<Self> {
        let expected_twist = self.curve().frobenius_twist_power(1).ok()?;
        if &expected_twist != codomain_twist {
            return None;
        }

        let rhs_factor = self.frobenius_quadratic_rhs_factor();
        let adjusted_b = self.b_part().div(&rhs_factor).ok()?;

        let a_part = self
            .a_part()
            .inverse_absolute_frobenius_pullback_from_twist()?;
        let b_part = adjusted_b.inverse_absolute_frobenius_pullback_from_twist()?;

        Some(Self::new(codomain_twist.clone(), a_part, b_part))
    }
}
