use num_complex::Complex64;

use super::super::AnalyticCurveError;
use crate::fields::ComplexApprox;
use crate::numerics::ApproxTolerance;

/// Three approximate roots of the cubic factor in
/// `4x^3 - g₂x - g₃ = 4(x - e_1)(x - e_2)(x - e_3)`.
///
/// This type stores one caller-supplied triple of roots and rejects triples
/// with approximately repeated entries under the requested tolerance.
///
/// The stored order is preserved from construction time, but that order does
/// not carry any canonical geometric or algebraic meaning.
#[derive(Clone, Debug, PartialEq)]
pub struct WeierstrassCubicRoots {
    roots: [Complex64; 3],
}

impl WeierstrassCubicRoots {
    /// Builds a validated triple of approximate cubic roots.
    ///
    /// The constructor preserves the caller-supplied order and rejects when
    /// any two roots are approximately equal under `tolerance`.
    pub fn new(
        e1: Complex64,
        e2: Complex64,
        e3: Complex64,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        let roots = Self {
            roots: [e1, e2, e3],
        };
        roots.validate_distinct(tolerance)?;
        Ok(roots)
    }

    /// Returns the stored triple of roots in construction order.
    pub fn roots(&self) -> [&Complex64; 3] {
        [&self.roots[0], &self.roots[1], &self.roots[2]]
    }

    /// Returns the first elementary symmetric sum `e₁ + e₂ + e₃`.
    pub fn sum(&self) -> Complex64 {
        self.roots[0] + self.roots[1] + self.roots[2]
    }

    /// Returns the second elementary symmetric sum
    /// `e₁e₂ + e₁e₃ + e₂e₃`.
    pub fn pairwise_products_sum(&self) -> Complex64 {
        self.roots[0] * self.roots[1]
            + self.roots[0] * self.roots[2]
            + self.roots[1] * self.roots[2]
    }

    /// Returns the product `e₁e₂e₃`.
    pub fn product(&self) -> Complex64 {
        self.roots[0] * self.roots[1] * self.roots[2]
    }

    /// Returns the coefficient of `x²` in
    /// `4(x - e₁)(x - e₂)(x - e₃)`.
    pub(crate) fn x_squared_coefficient(&self) -> Complex64 {
        -Complex64::new(4.0, 0.0) * self.sum()
    }

    /// Returns the implied invariant `g₂`.
    pub fn g2(&self) -> Complex64 {
        -Complex64::new(4.0, 0.0) * self.pairwise_products_sum()
    }

    /// Returns the implied invariant `g₃`.
    pub fn g3(&self) -> Complex64 {
        Complex64::new(4.0, 0.0) * self.product()
    }

    /// Returns the smallest pairwise distance among the three roots.
    pub fn min_pairwise_distance(&self) -> f64 {
        let d12 = (self.roots[0] - self.roots[1]).norm();
        let d13 = (self.roots[0] - self.roots[2]).norm();
        let d23 = (self.roots[1] - self.roots[2]).norm();

        d12.min(d13).min(d23)
    }

    /// Returns whether some pair of roots is approximately repeated under
    /// `tolerance`.
    pub fn has_repeated_root_approx(&self, tolerance: ApproxTolerance) -> bool {
        ComplexApprox::eq_with_tolerance(&self.roots[0], &self.roots[1], tolerance)
            || ComplexApprox::eq_with_tolerance(&self.roots[0], &self.roots[2], tolerance)
            || ComplexApprox::eq_with_tolerance(&self.roots[1], &self.roots[2], tolerance)
    }

    /// Returns one permutation witness sending this stored triple to
    /// `other`, when the two triples agree approximately up to reordering.
    #[cfg(test)]
    pub(crate) fn matching_permutation(
        &self,
        other: &Self,
        tolerance: ApproxTolerance,
    ) -> Option<[usize; 3]> {
        const PERMUTATIONS: [[usize; 3]; 6] = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];

        PERMUTATIONS.into_iter().find(|permutation| {
            ComplexApprox::eq_with_tolerance(
                &self.roots[0],
                &other.roots[permutation[0]],
                tolerance,
            ) && ComplexApprox::eq_with_tolerance(
                &self.roots[1],
                &other.roots[permutation[1]],
                tolerance,
            ) && ComplexApprox::eq_with_tolerance(
                &self.roots[2],
                &other.roots[permutation[2]],
                tolerance,
            )
        })
    }

    /// Returns whether the two root triples agree approximately up to
    /// permutation.
    #[cfg(test)]
    pub(crate) fn approximately_matches_up_to_permutation(
        &self,
        other: &Self,
        tolerance: ApproxTolerance,
    ) -> bool {
        self.matching_permutation(other, tolerance).is_some()
    }

    /// Returns whether the cubic is approximately depressed, equivalently
    /// whether `e₁ + e₂ + e₃ ≈ 0`.
    #[cfg(test)]
    pub(crate) fn is_approximately_depressed(&self, tolerance: ApproxTolerance) -> bool {
        ComplexApprox::is_zero_with_tolerance(&self.sum(), tolerance)
    }

    /// Validates that all three roots are pairwise distinct under `tolerance`.
    pub(crate) fn validate_distinct(
        &self,
        tolerance: ApproxTolerance,
    ) -> Result<(), AnalyticCurveError> {
        if self.has_repeated_root_approx(tolerance) {
            return Err(AnalyticCurveError::RepeatedCubicRoot);
        }

        Ok(())
    }
}
