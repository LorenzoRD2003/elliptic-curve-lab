use crate::elliptic_curves::analytic::AnalyticCurveError;

/// Truncation policy for first approximations to classical elliptic
/// functions attached to a lattice `Λ`.
///
/// The Weierstrass functions `℘`, `ζ`, and `σ` are defined by infinite
/// lattice expressions. This type provides a small validated knob for the
/// educational milestone where those expressions will first be approximated
/// by finite symmetric boxes of lattice indices.
///
/// Invariants:
/// - `radius > 0`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EllipticFunctionTruncation {
    radius: usize,
}

impl EllipticFunctionTruncation {
    /// Builds a validated truncation radius for elliptic-function
    /// approximations.
    ///
    /// Radius `0` is rejected because it would collapse the truncation to an
    /// origin-only box and stop being meaningful for lattice-based
    /// approximations.
    pub fn new(radius: usize) -> Result<Self, AnalyticCurveError> {
        if radius == 0 {
            return Err(AnalyticCurveError::InvalidTruncationRadius);
        }

        Ok(Self { radius })
    }

    /// Returns the stored truncation radius.
    pub fn radius(&self) -> usize {
        self.radius
    }

    /// Returns a small default radius intended for hand-checkable first
    /// experiments with truncated elliptic-function formulas.
    pub fn default_educational() -> Self {
        Self { radius: 2 }
    }
}
