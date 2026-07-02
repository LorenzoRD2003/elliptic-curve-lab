use core::cmp::Ordering;
use core::fmt;

use crate::elliptic_curves::frobenius::FrobeniusTrace;
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;
use num_bigint::{BigInt, BigUint};

/// Characteristic polynomial of the relative Frobenius `π_q`.
///
/// For an elliptic curve over `F_q` with Frobenius trace `t`, the relative
/// Frobenius satisfies `χ_{π_q}(T) = T^2 - tT + q`.
///
/// This value object stores exactly the finite base-field metadata and the
/// trace needed to write that polynomial explicitly.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusCharacteristicPolynomial {
    base_field: FiniteFieldDescriptor,
    trace: BigInt,
}

impl FrobeniusTrace {
    /// Returns the characteristic polynomial of the relative Frobenius `π_q`.
    ///
    /// If this trace package records `t = q + 1 - #E(F_q)`, then the returned
    /// polynomial is `χ_{π_q}(T) = T^2 - tT + q`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn characteristic_polynomial(&self) -> FrobeniusCharacteristicPolynomial {
        FrobeniusCharacteristicPolynomial::new(self.base_field().clone(), self.trace().clone())
    }
}

impl FrobeniusCharacteristicPolynomial {
    /// Builds `χ_{π_q}(T) = T^2 - tT + q` from `F_q` metadata and the trace `t`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn new(base_field: FiniteFieldDescriptor, trace: impl Into<BigInt>) -> Self {
        Self {
            base_field,
            trace: trace.into(),
        }
    }

    /// Returns the finite base-field descriptor for `F_q`.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        &self.base_field
    }

    /// Returns the finite base-field cardinality `q`.
    ///
    /// Complexity: `Θ(1)` after the descriptor-level `q = p^r` computation.
    pub fn field_order(&self) -> BigUint {
        self.base_field.cardinality_biguint()
    }

    /// Returns the Frobenius trace `t`.
    pub fn trace(&self) -> BigInt {
        self.trace.clone()
    }

    /// Returns the discriminant `t^2 - 4q`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn discriminant(&self) -> BigInt {
        &self.trace * &self.trace - BigInt::from(BigUint::from(4u8) * self.field_order())
    }

    /// Evaluates `χ_{π_q}(x) = x^2 - tx + q` at an integer `x`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn evaluate_at_integer(&self, x: i64) -> BigInt {
        let x = BigInt::from(x);
        &x * &x - &self.trace * &x + BigInt::from(self.field_order())
    }

    /// Returns a compact educational string such as `T^2 - 3T + 43`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn pretty(&self) -> String {
        let q = self.field_order();
        let linear_term = match self.trace.cmp(&BigInt::from(0u8)) {
            Ordering::Greater => format!(" - {}T", self.trace),
            Ordering::Less => format!(" + {}T", -&self.trace),
            Ordering::Equal => String::new(),
        };
        format!("T^2{linear_term} + {q}")
    }
}

impl fmt::Display for FrobeniusCharacteristicPolynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.pretty())
    }
}
