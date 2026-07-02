use crate::elliptic_curves::{
    CurveError,
    frobenius::{
        HasseInterval,
        torsion::{FrobeniusTorsionMatrixError, ModNMatrix2},
    },
    traits::RelativeFrobeniusCurveModel,
};
use crate::fields::traits::*;
use crate::fields::{finite_field_descriptor::FiniteFieldDescriptor, traits::FiniteField};
use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
use num_traits::{One, ToPrimitive, Zero};

/// Frobenius trace data recovered from a point count over a finite base field.
///
/// For an elliptic curve over `F_q`, the trace of the relative Frobenius
/// `π_q` is the integer `t = q + 1 - #E(F_q)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusTrace {
    /// Finite base-field metadata for `F_q`.
    base_field: FiniteFieldDescriptor,
    /// The counted order `#E(F_q)`.
    curve_order: BigUint,
    /// The Frobenius trace `t = q + 1 - #E(F_q)`.
    trace: BigInt,
}

impl FrobeniusTrace {
    /// Builds a validated Frobenius-trace package from `F_q` and `#E(F_q)`.
    ///
    /// Complexity: `Θ(1)`
    pub fn from_order(
        base_field: FiniteFieldDescriptor,
        curve_order: impl ToBigUint,
    ) -> Result<Self, CurveError> {
        let curve_order =
            curve_order
                .to_biguint()
                .ok_or_else(|| CurveError::InvalidCurveOrder {
                    order: BigUint::zero(),
                })?;
        let field_order = base_field.cardinality_biguint();

        if curve_order.is_zero() {
            return Err(CurveError::InvalidCurveOrder {
                order: BigUint::zero(),
            });
        }

        let trace = BigInt::from(field_order) + BigInt::one() - BigInt::from(curve_order.clone());

        Ok(Self {
            base_field,
            curve_order,
            trace,
        })
    }

    /// Reconstructs `#E(F_q)` from `F_q` and the Frobenius trace `t`.
    ///
    /// Complexity: `Θ(1)`
    pub fn curve_order_from_trace(
        base_field: FiniteFieldDescriptor,
        trace: impl ToBigInt,
    ) -> Result<BigUint, CurveError> {
        let trace = trace
            .to_bigint()
            .ok_or_else(|| CurveError::InvalidFrobeniusTrace {
                trace: BigInt::zero(),
            })?;
        let field_order = BigInt::from(base_field.cardinality_biguint());

        let curve_order = field_order + BigInt::one() - trace.clone();
        if curve_order <= BigInt::zero() {
            return Err(CurveError::InvalidFrobeniusTrace { trace });
        }

        curve_order
            .to_biguint()
            .ok_or(CurveError::InvalidFrobeniusTrace { trace })
    }

    /// Returns the finite base-field descriptor for `F_q`.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        &self.base_field
    }

    /// Returns the finite base-field cardinality `q`.
    pub fn field_order(&self) -> BigUint {
        self.base_field.cardinality_biguint()
    }

    pub fn curve_order(&self) -> BigUint {
        self.curve_order.clone()
    }

    pub fn trace(&self) -> BigInt {
        self.trace.clone()
    }

    /// Returns the discrete Hasse interval `H(q)` attached to this base field.
    ///
    /// If this trace package stores data for a curve over `F_q`, the returned
    /// interval is
    ///
    /// `H(q) = [ceil(q + 1 - 2 sqrt(q)), floor(q + 1 + 2 sqrt(q))]`,
    ///
    /// the standard integer search interval that must contain `#E(F_q)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn hasse_interval(&self) -> HasseInterval {
        HasseInterval::from_trace(self)
    }

    pub(crate) fn assert_compatible_with_curve<E>(&self) -> Result<(), CurveError>
    where
        E: RelativeFrobeniusCurveModel + ?Sized,
        E::BaseField: FiniteField,
    {
        Self::assert_base_field_compatible_with_curve::<E>(self.base_field())
    }

    pub(crate) fn assert_base_field_compatible_with_curve<E>(
        base_field: &FiniteFieldDescriptor,
    ) -> Result<(), CurveError>
    where
        E: RelativeFrobeniusCurveModel + ?Sized,
        E::BaseField: FiniteField,
    {
        let characteristic = E::BaseField::characteristic().to_biguint();
        let curve_base_field =
            FiniteFieldDescriptor::new(characteristic, E::BaseField::extension_degree()).expect(
                "finite field implementations should expose internally consistent metadata",
            );
        if &curve_base_field != base_field {
            return Err(CurveError::IncompatibleFrobeniusBaseField {
                curve_characteristic: curve_base_field.characteristic.clone(),
                curve_extension_degree: curve_base_field.extension_degree.get(),
                polynomial_characteristic: base_field.characteristic.clone(),
                polynomial_extension_degree: base_field.extension_degree.get(),
            });
        }
        Ok(())
    }

    /// Checks whether the trace of a Frobenius matrix on `E[n]` agrees with
    /// this stored Frobenius trace modulo `n`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn trace_matches_torsion_matrix_mod_n(
        &self,
        matrix: &ModNMatrix2,
    ) -> Result<bool, FrobeniusTorsionMatrixError> {
        let modulus = BigInt::from(matrix.modulus());
        let trace = self.trace();
        let reduced_trace = ((trace % &modulus) + &modulus) % &modulus;
        let matrix_trace = BigInt::from(matrix.trace_mod_n());
        Ok(reduced_trace == matrix_trace)
    }

    /// Checks whether the determinant of a Frobenius matrix on `E[n]` agrees
    /// with `q mod n`, where this trace package is for a curve over `F_q`.
    ///
    /// Complexity: `Θ(1)` plus the determinant reduction of `matrix`.
    pub fn determinant_matches_torsion_matrix_mod_n(
        &self,
        matrix: &ModNMatrix2,
    ) -> Result<bool, FrobeniusTorsionMatrixError> {
        let determinant = matrix.determinant_mod_n()?;
        let reduced_field_order = (&self.field_order() % BigUint::from(matrix.modulus()))
            .to_usize()
            .expect("reduction modulo usize should fit in usize");
        Ok(reduced_field_order == determinant)
    }
}

pub(crate) fn curve_order_from_field_order_and_trace(
    field_order: &BigUint,
    trace: &BigInt,
) -> Result<BigUint, CurveError> {
    let curve_order = BigInt::from(field_order.clone()) + BigInt::one() - trace;
    curve_order
        .to_biguint()
        .ok_or_else(|| CurveError::InvalidHasseIntervalFieldOrder {
            field_order: field_order.clone(),
        })
}
