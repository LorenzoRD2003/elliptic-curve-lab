use num_bigint::BigInt;
use num_complex::Complex64;
use num_rational::BigRational;
use num_traits::ToPrimitive;

use crate::elliptic_curves::analytic::{AnalyticCurveError, q_expansion::QExpansionTruncation};

/// Stored coefficient table for a modular `q`-expansion.
///
/// This value object stores a finite exact prefix
/// `a_0 q^m + a_1 q^(m+1) + ... + a_r q^(m+r)`
/// through the starting exponent `m`, and an exact rational
/// coefficient list `[a_0, ..., a_r]`
///
/// It is intentionally neutral about whether those coefficients happen to be
/// integral. That matters for the current analytic stages because:
///
/// - the shipped `j(q)` table has integer coefficients
/// - the holomorphic Eisenstein family `E_k(q)` has exact rational
///   coefficients in general, for example at weight `12`
///
/// So this type is the common exact coefficient layer shared by both modular
/// functions and modular forms before any conversion to approximate
/// `Complex64` evaluation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ModularQExpansionCoefficients {
    start_exponent: i32,
    coefficients: Vec<BigRational>,
}

impl ModularQExpansionCoefficients {
    /// Builds a coefficient table from an initial exponent and explicit
    /// exact rational coefficients.
    ///
    /// The stored coefficients correspond to
    /// `q^start_exponent, q^(start_exponent+1), ...`.
    pub(crate) fn new(start_exponent: i32, coefficients: Vec<BigRational>) -> Self {
        Self {
            start_exponent,
            coefficients,
        }
    }

    /// Builds a coefficient table from exact integer coefficients.
    ///
    /// This is a convenience constructor for families such as the current
    /// educational `j(q)` table whose stored prefix happens to lie in `ℤ`.
    pub(crate) fn from_integers<I, C>(start_exponent: i32, coefficients: I) -> Self
    where
        I: IntoIterator<Item = C>,
        C: Into<BigInt>,
    {
        Self::new(
            start_exponent,
            coefficients
                .into_iter()
                .map(|coefficient| BigRational::from_integer(coefficient.into()))
                .collect(),
        )
    }

    /// Returns the first exponent represented by this table.
    #[cfg(test)]
    pub(crate) fn start_exponent(&self) -> i32 {
        self.start_exponent
    }

    /// Returns the last exponent represented by this table, if any
    /// coefficients are present.
    #[cfg(test)]
    pub(crate) fn end_exponent(&self) -> Option<i32> {
        self.coefficients
            .len()
            .checked_sub(1)
            .map(|offset| self.start_exponent + offset as i32)
    }

    /// Returns the stored exact rational coefficients.
    #[cfg(test)]
    pub(crate) fn coefficients(&self) -> &[BigRational] {
        &self.coefficients
    }

    /// Returns the number of stored coefficients.
    pub(crate) fn len(&self) -> usize {
        self.coefficients.len()
    }

    /// Returns whether the table is empty.
    #[cfg(test)]
    pub(crate) fn is_empty(&self) -> bool {
        self.coefficients.is_empty()
    }

    /// Returns the coefficient of `q^exponent`, if this table stores it.
    #[cfg(test)]
    pub(crate) fn coefficient_of(&self, exponent: i32) -> Option<BigRational> {
        let offset = exponent.checked_sub(self.start_exponent)?;
        let index = usize::try_from(offset).ok()?;
        self.coefficients.get(index).cloned()
    }

    /// Returns the shorter exact prefix determined by the supplied truncation.
    pub(crate) fn truncated(
        &self,
        truncation: QExpansionTruncation,
    ) -> Result<ModularQExpansionCoefficients, AnalyticCurveError> {
        if truncation.terms() > self.len() {
            return Err(AnalyticCurveError::InvalidSeriesPrecision);
        }

        Ok(Self::new(
            self.start_exponent,
            self.coefficients
                .iter()
                .take(truncation.terms())
                .cloned()
                .collect(),
        ))
    }

    /// Evaluates the full stored finite sum `Σ a_k q^(start_exponent + k)`.
    ///
    /// The stored rational table remains exact; the conversion to floating-point
    /// complex arithmetic happens only at this evaluation boundary.
    ///
    /// TODO:
    /// if later stages need algebraic operations on `q`-series themselves,
    /// promote this value object toward a small formal-series layer instead of
    /// using only direct numerical evaluation into `Complex64`.
    pub(crate) fn evaluate_at(&self, q: Complex64) -> Result<Complex64, AnalyticCurveError> {
        self.coefficients.iter().enumerate().try_fold(
            Complex64::new(0.0, 0.0),
            |accumulator, (offset, coefficient)| {
                let exponent = self.start_exponent + offset as i32;
                let coefficient_as_f64 = coefficient
                    .to_f64()
                    .ok_or(AnalyticCurveError::NumericalComparisonFailed)?;

                Ok(accumulator
                    + Complex64::new(coefficient_as_f64, 0.0) * pow_complex_i32(q, exponent))
            },
        )
    }

    /// Evaluates the truncated partial sum
    /// `Σ_{k=0}^{N-1} a_k q^(start_exponent + k)`
    /// at the supplied modular parameter `q`, where `N = truncation.terms()`.
    #[cfg(test)]
    pub(crate) fn evaluate_truncated_at(
        &self,
        q: Complex64,
        truncation: QExpansionTruncation,
    ) -> Result<Complex64, AnalyticCurveError> {
        self.truncated(truncation)?.evaluate_at(q)
    }
}

fn pow_complex_i32(base: Complex64, exponent: i32) -> Complex64 {
    if exponent >= 0 {
        base.powu(exponent as u32)
    } else {
        Complex64::new(1.0, 0.0) / base.powu(exponent.unsigned_abs())
    }
}
