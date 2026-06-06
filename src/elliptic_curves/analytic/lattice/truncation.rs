use crate::elliptic_curves::analytic::AnalyticCurveError;

/// Truncation policy for finite square-box approximations to sums over a
/// complex lattice `Λ = ℤω₁ + ℤω₂`.
///
/// Many analytic constructions attached to elliptic curves involve infinite
/// sums over all lattice points `mω₁ + nω₂`. This type models the simple
/// truncation that keeps only the integer index box `-r ≤ m ≤ r`, `-r ≤ n ≤ r`.
///
/// Some routines, such as classical Eisenstein or Weierstrass lattice
/// sums, use the punctured box obtained by removing the origin. This type
/// therefore exposes both the full box count and the nonzero box count
/// explicitly.
///
/// Invariants:
/// - `radius > 0`
/// - the corresponding square-box term count `(2r + 1)²` fits in `usize`
///
/// Note: This is **not** a geometric disc or norm-bound truncation in `ℂ`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LatticeSumTruncation {
    radius: usize,
}

impl LatticeSumTruncation {
    /// Builds a validated square-box truncation radius.
    ///
    /// The radius must be positive. Radius `0` is rejected because the
    /// resulting box would contain only the origin.
    ///
    /// This constructor also rejects radii whose square-box term count
    /// `(2r + 1)²` would overflow `usize`.
    pub fn new(radius: usize) -> Result<Self, AnalyticCurveError> {
        if radius == 0 || Self::terms_in_square_box_for_radius(radius).is_none() {
            return Err(AnalyticCurveError::InvalidTruncationRadius);
        }

        Ok(Self { radius })
    }

    /// Returns the stored truncation radius `r`.
    pub fn radius(&self) -> usize {
        self.radius
    }

    /// Returns the number of index pairs in the square box
    /// `-r ≤ m ≤ r`, `-r ≤ n ≤ r`.
    ///
    /// Equivalently, this is `(2r + 1)²`, including the origin `(0, 0)`.
    pub fn terms_in_square_box(&self) -> usize {
        Self::terms_in_square_box_for_radius(self.radius)
            .expect("validated truncation radius must have a finite square-box term count")
    }

    /// Returns the number of nonzero index pairs in the square box
    /// `-r ≤ m ≤ r`, `-r ≤ n ≤ r`.
    ///
    /// Equivalently, this is `(2r + 1)² - 1`.
    pub fn nonzero_terms_in_square_box(&self) -> usize {
        self.terms_in_square_box()
            .checked_sub(1)
            .expect("validated positive radius box must contain the origin")
    }

    /// Returns a small default truncation intended for hand-checkable
    /// examples and first experiments.
    pub fn default_educational() -> Self {
        Self { radius: 2 }
    }

    /// Returns a somewhat larger truncation intended for side-by-side
    /// comparison against [`Self::default_educational()`].
    pub fn larger_for_comparison() -> Self {
        Self { radius: 4 }
    }

    fn terms_in_square_box_for_radius(radius: usize) -> Option<usize> {
        let side_length = radius.checked_mul(2)?.checked_add(1)?;
        side_length.checked_mul(side_length)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::analytic::lattice::LatticeSumTruncation;
    use crate::elliptic_curves::analytic::AnalyticCurveError;

    #[test]
    fn truncation_requires_positive_radius() {
        assert_eq!(
            LatticeSumTruncation::new(0),
            Err(AnalyticCurveError::InvalidTruncationRadius)
        );
    }

    #[test]
    fn truncation_exposes_radius_and_square_box_count() {
        let truncation = LatticeSumTruncation::new(3).unwrap();

        assert_eq!(truncation.radius(), 3);
        assert_eq!(truncation.terms_in_square_box(), 49);
        assert_eq!(truncation.nonzero_terms_in_square_box(), 48);
    }

    #[test]
    fn educational_presets_are_explicit_and_ordered() {
        let educational = LatticeSumTruncation::default_educational();
        let larger = LatticeSumTruncation::larger_for_comparison();

        assert_eq!(educational.radius(), 2);
        assert_eq!(educational.terms_in_square_box(), 25);
        assert_eq!(educational.nonzero_terms_in_square_box(), 24);
        assert_eq!(larger.radius(), 4);
        assert_eq!(larger.terms_in_square_box(), 81);
        assert_eq!(larger.nonzero_terms_in_square_box(), 80);
        assert!(larger.radius() > educational.radius());
    }

    #[test]
    fn truncation_rejects_radius_when_square_box_count_would_overflow() {
        let overflowing_radius = (usize::MAX / 2) + 1;

        assert_eq!(
            LatticeSumTruncation::new(overflowing_radius),
            Err(AnalyticCurveError::InvalidTruncationRadius)
        );
    }
}
