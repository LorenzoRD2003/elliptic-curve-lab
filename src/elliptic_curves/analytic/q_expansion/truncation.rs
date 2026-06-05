use crate::elliptic_curves::analytic::AnalyticCurveError;

/// Truncation policy for short educational `q`-expansions.
///
/// The stored `terms` counts how many coefficients from the chosen family are
/// retained in order, starting from that family's first exponent.
///
/// Examples:
///
/// - for `E₄(q) = 1 + 240q + 2160q² + ...`, `terms = 1` means only the
///   constant term `1`
/// - for `j(q) = q^{-1} + 744 + 196884q + ...`, `terms = 2` means `q^{-1} + 744`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QExpansionTruncation {
    terms: usize,
}

impl QExpansionTruncation {
    /// Builds a validated `q`-expansion truncation.
    ///
    /// The current educational implementation requires `terms >= 1`.
    pub fn new(terms: usize) -> Result<Self, AnalyticCurveError> {
        if terms >= 1 {
            Ok(Self { terms })
        } else {
            Err(AnalyticCurveError::InvalidSeriesPrecision)
        }
    }

    /// Returns the number of stored nonnegative-power terms.
    pub fn terms(&self) -> usize {
        self.terms
    }

    /// Returns a compact default truncation suitable for first experiments.
    pub fn default_educational() -> Self {
        Self { terms: 4 }
    }
}
