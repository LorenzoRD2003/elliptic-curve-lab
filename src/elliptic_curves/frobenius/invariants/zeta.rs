use core::cmp::Ordering;
use core::fmt;

use crate::elliptic_curves::frobenius::{FrobeniusCharacteristicPolynomial, FrobeniusTrace};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;
use num_bigint::{BigInt, BigUint};

/// Local zeta function of an elliptic curve over a finite field `F_q`.
///
/// If the relative Frobenius `π_q` has trace `t`, then the local zeta
/// function is
///
/// `Z(E/F_q, T) = (1 - tT + qT^2) / ((1 - T)(1 - qT))`.
///
/// In the current educational implementation, this object is derived from the
/// characteristic polynomial `χ_{π_q}(T) = T^2 - tT + q`, so the zeta
/// numerator and the Frobenius characteristic polynomial stay tied to one
/// source of truth.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusLocalZetaFunction {
    characteristic_polynomial: FrobeniusCharacteristicPolynomial,
}

impl FrobeniusTrace {
    /// Returns the local zeta function `Z(E/F_q, T)` attached to this trace.
    ///
    /// If this trace package records `t = q + 1 - #E(F_q)`, then the returned
    /// zeta function is
    /// `Z(E/F_q, T) = (1 - tT + qT^2) / ((1 - T)(1 - qT))`.
    pub fn local_zeta_function(&self) -> FrobeniusLocalZetaFunction {
        self.characteristic_polynomial().local_zeta_function()
    }
}

impl FrobeniusCharacteristicPolynomial {
    /// Returns the local zeta function determined by `χ_{π_q}(T)`.
    ///
    /// If `χ_{π_q}(T) = T^2 - tT + q`, then the returned zeta function is
    /// `Z(E/F_q, T) = (1 - tT + qT^2) / ((1 - T)(1 - qT))`.
    pub fn local_zeta_function(&self) -> FrobeniusLocalZetaFunction {
        FrobeniusLocalZetaFunction::from_characteristic_polynomial(self.clone())
    }
}

impl FrobeniusLocalZetaFunction {
    /// Builds the local zeta function from `χ_{π_q}(T) = T^2 - tT + q`.
    pub fn from_characteristic_polynomial(
        characteristic_polynomial: FrobeniusCharacteristicPolynomial,
    ) -> Self {
        Self {
            characteristic_polynomial,
        }
    }

    /// Returns the underlying characteristic polynomial `χ_{π_q}(T)`.
    pub fn characteristic_polynomial(&self) -> &FrobeniusCharacteristicPolynomial {
        &self.characteristic_polynomial
    }

    /// Returns the finite base-field descriptor for `F_q`.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        self.characteristic_polynomial.base_field()
    }

    /// Returns the finite base-field cardinality `q`.
    pub fn field_order(&self) -> BigUint {
        self.characteristic_polynomial.field_order()
    }

    /// Returns the Frobenius trace `t`.
    pub fn trace(&self) -> BigInt {
        self.characteristic_polynomial.trace()
    }

    /// Returns a compact educational string for the local zeta function.
    ///
    /// Example: `Z(E/F_17, T) = (1 - 2T + 17T²) / ((1 - T)(1 - 17T))`.
    pub fn pretty(&self) -> String {
        format!(
            "Z(E/{}, T) = ({}) / ({})",
            self.base_field(),
            self.numerator_pretty(),
            self.denominator_pretty()
        )
    }

    /// Returns the numerator `1 - tT + qT²` in a compact educational form.
    pub(crate) fn numerator_pretty(&self) -> String {
        let trace = self.trace();
        let field_order = self.field_order();
        let linear_term = match trace.cmp(&BigInt::from(0u8)) {
            Ordering::Greater => format!(" - {trace}T"),
            Ordering::Less => format!(" + {}T", -trace),
            Ordering::Equal => String::new(),
        };

        format!("1{linear_term} + {field_order}T²")
    }

    /// Returns the denominator `(1 - T)(1 - qT)` in a compact educational form.
    pub(crate) fn denominator_pretty(&self) -> String {
        format!("(1 - T)(1 - {}T)", self.field_order())
    }
}

impl fmt::Display for FrobeniusLocalZetaFunction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.pretty())
    }
}
