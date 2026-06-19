use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::schoof::ReducedCurveQuotient,
};
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

use super::ReducedEndomorphism;

/// Result of one computation in the additive arithmetic of reduced
/// endomorphisms.
///
/// This keeps the additive zero endomorphism explicit, because the constant
/// map `P ↦ O` is not representable by one affine reduced pair
/// `(a(x), b(x) y)`.
#[derive(Debug)]
pub enum ReducedEndomorphismAdditiveResult<F: FiniteField> {
    Zero,
    Value(ReducedEndomorphism<F>),
    NonUnitDenominator { witness_gcd: DensePolynomial<F> },
}

impl<F: FiniteField> Clone for ReducedEndomorphismAdditiveResult<F> {
    fn clone(&self) -> Self {
        match self {
            Self::Zero => Self::Zero,
            Self::Value(value) => Self::Value(value.clone()),
            Self::NonUnitDenominator { witness_gcd } => Self::NonUnitDenominator {
                witness_gcd: witness_gcd.clone(),
            },
        }
    }
}

impl<F: FiniteField> PartialEq for ReducedEndomorphismAdditiveResult<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Zero, Self::Zero) => true,
            (Self::Value(lhs), Self::Value(rhs)) => lhs == rhs,
            (
                Self::NonUnitDenominator {
                    witness_gcd: lhs_witness,
                },
                Self::NonUnitDenominator {
                    witness_gcd: rhs_witness,
                },
            ) => lhs_witness == rhs_witness,
            _ => false,
        }
    }
}

impl<F: FiniteField> ReducedEndomorphismAdditiveResult<F> {
    /// Negates one additive result in the reduced endomorphism arithmetic.
    ///
    /// The additive zero and non-unit-denominator branches stay unchanged,
    /// while an affine representative `(a(x), b(x) y)` is replaced by its
    /// additive inverse `(a(x), -b(x) y)`.
    pub(crate) fn additive_inverse(self) -> Self {
        match self {
            Self::Zero => Self::Zero,
            Self::Value(value) => Self::Value(value.additive_inverse()),
            Self::NonUnitDenominator { witness_gcd } => Self::NonUnitDenominator { witness_gcd },
        }
    }

    /// Adds one non-zero affine reduced endomorphism to this additive result.
    ///
    /// This helper keeps the additive zero branch and the non-unit-denominator
    /// branch explicit, while delegating the genuine affine-affine addition to
    /// the short-Weierstrass group-law helper.
    pub(crate) fn combine_with_value(
        self,
        curve: &ShortWeierstrassCurve<F>,
        quotient: &ReducedCurveQuotient<F>,
        value: ReducedEndomorphism<F>,
    ) -> Result<Self, DensePolynomial<F>> {
        match self {
            Self::Zero => Ok(Self::Value(value)),
            Self::Value(accumulator_value) => curve
                .add_reduced_endomorphisms(quotient, &accumulator_value, &value)
                .into_result(),
            Self::NonUnitDenominator { witness_gcd } => Err(witness_gcd),
        }
    }

    /// Appends another additive result to this one.
    pub(crate) fn combine(
        self,
        curve: &ShortWeierstrassCurve<F>,
        quotient: &ReducedCurveQuotient<F>,
        rhs: Self,
    ) -> Result<Self, DensePolynomial<F>> {
        match rhs {
            Self::Zero => Ok(self),
            Self::Value(value) => self.combine_with_value(curve, quotient, value),
            Self::NonUnitDenominator { witness_gcd } => Err(witness_gcd),
        }
    }

    pub(crate) fn into_result(self) -> Result<Self, DensePolynomial<F>> {
        match self {
            Self::NonUnitDenominator { witness_gcd } => Err(witness_gcd),
            other => Ok(other),
        }
    }
}
