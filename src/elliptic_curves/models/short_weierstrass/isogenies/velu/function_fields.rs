use std::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    short_weierstrass::function_fields::{
        ShortWeierstrassFunction, ShortWeierstrassFunctionField, ShortWeierstrassFunctionFieldPoint,
    },
    short_weierstrass::isogenies::{
        VeluIsogeny,
        function_field_maps::{DifferentialPullbackReport, ShortWeierstrassFunctionFieldMap},
    },
};
use crate::fields::traits::Field;
use crate::isogenies::error::IsogenyError;

impl<F: Field + Clone> VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    /// Returns the pullback `φ*(x')` of the codomain `x`-coordinate.
    ///
    /// For the current short-Weierstrass Vélu normalization and `G* = G \ {O}`,
    /// this is the function `x + Σ_{Q ∈ G*} (x(P + Q) - x(Q))` on a generic
    /// point `P = (x, y)` of the domain curve.
    pub fn x_pullback(&self) -> ShortWeierstrassFunction<F> {
        let field = ShortWeierstrassFunctionField::<F>::new(self.domain.clone());
        let x = field.x();

        self.kernel_nonzero_points()
            .iter()
            .fold(x.clone(), |accumulator, kernel_point| {
                let translated_point = self
                    .translated_generic_point(kernel_point)
                    .expect("generic translation by a non-zero kernel point should be defined in the function field");
                let translated_x = translated_point
                    .x()
                    .expect("translation by a finite point stays affine");
                let kernel_x = ShortWeierstrassFunction::<F>::from_finite_affine_point(
                    self.domain.clone(),
                    kernel_point,
                    true,
                );
                let correction = translated_x
                    .sub(&kernel_x)
                    .expect("same-curve subtraction should work");

                accumulator
                    .add(&correction)
                    .expect("same-curve addition should work")
            })
    }

    /// Returns the pullback `φ*(y')` of the codomain `y`-coordinate.
    ///
    /// For the current short-Weierstrass Vélu normalization and `G* = G \ {O}`,
    /// this is the function `y + Σ_{Q ∈ G*} (y(P + Q) - y(Q))` on a generic
    /// point `P = (x, y)` of the domain curve.
    pub fn y_pullback(&self) -> ShortWeierstrassFunction<F> {
        let field = ShortWeierstrassFunctionField::<F>::new(self.domain.clone());
        let y = field.y();

        self.kernel_nonzero_points()
            .iter()
            .fold(y.clone(), |accumulator, kernel_point| {
                let translated_point = self
                    .translated_generic_point(kernel_point)
                    .expect("generic translation by a non-zero kernel point should be defined in the function field");
                let translated_y = translated_point
                    .y()
                    .expect("translation by a finite point stays affine");
                let kernel_y = ShortWeierstrassFunction::<F>::from_finite_affine_point(
                    self.domain.clone(),
                    kernel_point,
                    false,
                );
                let correction = translated_y
                    .sub(&kernel_y)
                    .expect("same-curve subtraction should work");

                accumulator
                    .add(&correction)
                    .expect("same-curve addition should work")
            })
    }

    /// Returns the current Vélu map as a pullback `φ* : F(E') -> F(E)`.
    pub fn as_function_field_map(&self) -> ShortWeierstrassFunctionFieldMap<F> {
        ShortWeierstrassFunctionFieldMap::new(
            self.domain.clone(),
            self.codomain.clone(),
            self.x_pullback(),
            self.y_pullback(),
        )
        .expect("Vélu pullbacks should satisfy the codomain equation")
    }

    /// Returns the current differential pullback report for this Vélu isogeny.
    pub fn differential_pullback_report(
        &self,
    ) -> Result<DifferentialPullbackReport<F>, IsogenyError> {
        self.as_function_field_map().differential_pullback_report()
    }

    fn translated_generic_point(
        &self,
        kernel_point: &AffinePoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        let field = ShortWeierstrassFunctionField::<F>::new(self.domain.clone());
        field.translate_generic_point_by_base_point(kernel_point)
    }
}
