use std::hash::Hash;

use crate::elliptic_curves::{
    frobenius::torsion::matrix::{FrobeniusTorsionMatrixError, TorsionCoordinateCurveModel},
    traits::FiniteGroupCurveModel,
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NTorsionBasis<P> {
    n: usize,
    first: P,
    second: P,
}

impl<P> NTorsionBasis<P> {
    pub(crate) fn unchecked(n: usize, first: P, second: P) -> Self {
        Self { n, first, second }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn first(&self) -> &P {
        &self.first
    }

    pub fn second(&self) -> &P {
        &self.second
    }
}

impl<P: Clone + Eq + Hash> NTorsionBasis<P> {
    pub fn new<C: FiniteGroupCurveModel<Point = P>>(
        curve: &C,
        n: usize,
        first: P,
        second: P,
    ) -> Result<Self, FrobeniusTorsionMatrixError>
    where
        C::BaseField: EnumerableFiniteField + SqrtField,
    {
        curve.validate_basis(n, &first, &second)?;
        Ok(Self { n, first, second })
    }
}
