use std::hash::Hash;

use crate::elliptic_curves::short_weierstrass::{
    ShortWeierstrassCurve, function_fields::ShortWeierstrassFunctionField,
    isogenies::function_field_maps::ShortWeierstrassFunctionFieldMap,
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{error::IsogenyError, scalar_multiplication::ScalarMultiplicationIsogeny};

impl<F> ScalarMultiplicationIsogeny<ShortWeierstrassCurve<F>>
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    /// Returns the pullback map `[n]^* : F(E) -> F(E)` induced by scalar
    /// multiplication on the generic point.
    ///
    /// Let `P_gen = (x, y)` be the generic point of `E` viewed as a point of
    /// `E(F(E))`. Then the multiplication-by-`n` map is determined by the
    /// coordinates of `[n]P_gen = (X_n, Y_n)`.
    ///
    /// This method computes that generic multiple inside the existing
    /// short-Weierstrass function-field layer and returns the pullback
    /// `[n]^*(x) = X_n`, `[n]^*(y) = Y_n`.
    ///
    /// Since the constructor of [`ScalarMultiplicationIsogeny`] rejects the
    /// zero scalar, the image of the generic point is expected to stay affine
    /// in the current short-Weierstrass presentation.
    ///
    /// Complexity: `Θ(log n)` generic-point additions/doublings in `E(F(E))`
    /// from the double-and-add ladder, plus one final pullback-map
    /// validation.
    pub fn as_function_field_map(
        &self,
    ) -> Result<ShortWeierstrassFunctionFieldMap<F>, IsogenyError> {
        let field = ShortWeierstrassFunctionField::<F>::new(self.curve().clone());
        let image = field.generic_point_multiple(self.scalar())?;

        ShortWeierstrassFunctionFieldMap::new(
            self.curve().clone(),
            self.curve().clone(),
            image.x().unwrap().clone(),
            image.y().unwrap().clone(),
        )
    }
}
