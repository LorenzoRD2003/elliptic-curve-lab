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

/// Educational report for ordinary versus supersingular classification.
///
/// This report is intentionally derived from [`FrobeniusTrace`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusCurveTypeReport {
    frobenius_trace: FrobeniusTrace,
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
        let characteristic = i128::from(self.base_field().characteristic);
        let trace = i128::from(self.trace());
        if trace.rem_euclid(characteristic) == 0 {
            FrobeniusCurveType::Supersingular
        } else {
            FrobeniusCurveType::Ordinary
        }
    }

    /// Returns the ordinary/supersingular classification report derived from this trace.
    pub fn curve_type_report(&self) -> FrobeniusCurveTypeReport {
        FrobeniusCurveTypeReport {
            frobenius_trace: self.clone(),
        }
    }
}

impl FrobeniusCurveTypeReport {
    /// Returns the Frobenius trace package used for the classification.
    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    /// Returns the ordinary/supersingular classification.
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
    pub fn trace_mod_characteristic(&self) -> u64 {
        let characteristic = self.frobenius_trace.base_field().characteristic;
        let trace = i128::from(self.frobenius_trace.trace());
        let characteristic_i128 = i128::from(characteristic);
        let residue = trace.rem_euclid(characteristic_i128);
        u64::try_from(residue).expect("trace mod characteristic should fit into u64")
    }

    /// Returns whether the base-field characteristic divides the Frobenius trace.
    pub fn characteristic_divides_trace(&self) -> bool {
        self.trace_mod_characteristic() == 0
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
