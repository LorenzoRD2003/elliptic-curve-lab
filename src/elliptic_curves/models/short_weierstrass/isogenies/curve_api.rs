use crate::fields::traits::*;
use std::collections::HashSet;
use std::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::isogenies::{
        VeluIsogeny,
        frobenius::{AbsoluteFrobeniusIsogeny, RelativeFrobeniusIsogeny},
    },
};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::{error::IsogenyError, scalar_multiplication::ScalarMultiplicationIsogeny};

impl<F: Field + Clone> ShortWeierstrassCurve<F>
where
    F::Elem: Clone + Eq + Hash,
{
    /// Builds a short-Weierstrass Vélu isogeny from an explicit finite kernel.
    pub fn velu_isogeny_from_points(
        &self,
        kernel_points: HashSet<AffinePoint<F>>,
    ) -> Result<VeluIsogeny<Self>, IsogenyError> {
        VeluIsogeny::from_points(self.clone(), kernel_points)
    }
}

impl<F: EnumerableFiniteField + SqrtField + Clone> ShortWeierstrassCurve<F>
where
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    /// Builds a short-Weierstrass Vélu isogeny from one cyclic kernel generator.
    pub fn velu_isogeny_from_generator(
        &self,
        generator: AffinePoint<F>,
    ) -> Result<VeluIsogeny<Self>, IsogenyError> {
        VeluIsogeny::from_generator(self.clone(), generator)
    }

    /// Builds the scalar-multiplication isogeny `[n] : E -> E`.
    pub fn scalar_multiplication_isogeny(
        &self,
        scalar: u64,
    ) -> Result<ScalarMultiplicationIsogeny<Self>, IsogenyError> {
        ScalarMultiplicationIsogeny::new(self.clone(), scalar)
    }
}

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Builds the absolute Frobenius isogeny `Frob_p: E -> E^(p)`.
    pub fn absolute_frobenius_isogeny(&self) -> Result<AbsoluteFrobeniusIsogeny<F>, IsogenyError> {
        AbsoluteFrobeniusIsogeny::new(self.clone())
    }

    /// Builds the relative Frobenius isogeny `Frob_q: E -> E`.
    pub fn relative_frobenius_isogeny(&self) -> Result<RelativeFrobeniusIsogeny<F>, IsogenyError> {
        RelativeFrobeniusIsogeny::new(self.clone())
    }
}
