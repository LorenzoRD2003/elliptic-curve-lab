use num_bigint::BigInt;
use num_bigint::BigUint;
use num_traits::Zero;

use crate::elliptic_curves::frobenius::FrobeniusTrace;

/// Frobenius-side classification of an elliptic curve over a finite field.
///
/// For an elliptic curve over `F_q`, let `t` be the trace of the relative
/// Frobenius `π_q`, and let `p` be the characteristic of the base field.
///
/// - the curve is `Supersingular` if `p | t`
/// - the curve is `Ordinary` otherwise
///
/// In the prime-field case `F_p` with `p >= 5`, Hasse's bound forces
/// `|t| < p`, so `p | t` is equivalent to `t = 0`. We keep the divisibility
/// criterion as the primary API because it remains correct over extensions
/// `F_{p^n}` as well.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrobeniusCurveType {
    Ordinary,
    Supersingular,
}

impl FrobeniusTrace {
    /// Classifies the curve as ordinary or supersingular from the Frobenius trace.
    ///
    /// If the trace is `t` and the base-field characteristic is `p`, the
    /// classification criterion is:
    ///
    /// - `Supersingular` iff `p | t`
    /// - `Ordinary` iff `p ∤ t`
    ///
    /// Complexity: `Θ(1)`.
    pub fn curve_type(&self) -> FrobeniusCurveType {
        if self.characteristic_divides_trace() {
            FrobeniusCurveType::Supersingular
        } else {
            FrobeniusCurveType::Ordinary
        }
    }

    /// Returns the canonical residue class of `t mod p` in `{0, ..., p - 1}`.
    ///
    /// Here `t` is the Frobenius trace and `p` is the prime characteristic of
    /// the base field.
    pub fn trace_mod_characteristic(&self) -> BigUint {
        let characteristic = self.base_field().characteristic.clone();
        let characteristic_bigint = BigInt::from(characteristic.clone());
        let residue = (self.trace() % &characteristic_bigint + &characteristic_bigint)
            % &characteristic_bigint;
        residue
            .to_biguint()
            .expect("least nonnegative residue modulo p is nonnegative")
    }

    /// Returns whether the base-field characteristic divides the Frobenius trace.
    pub fn characteristic_divides_trace(&self) -> bool {
        self.trace_mod_characteristic().is_zero()
    }

    /// Returns whether the curve is ordinary.
    pub fn is_ordinary(&self) -> bool {
        self.curve_type() == FrobeniusCurveType::Ordinary
    }

    /// Returns whether the curve is supersingular.
    pub fn is_supersingular(&self) -> bool {
        self.curve_type() == FrobeniusCurveType::Supersingular
    }
}
