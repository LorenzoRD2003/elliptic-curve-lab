use std::hash::Hash;

use crate::elliptic_curves::traits::CurveModel;
use crate::isogenies::kernel::IsogenyKernel;

/// Educational scaffold for Vélu-style isogenies.
///
/// This type is represents an explicit isogeny constructed from a finite subgroup
/// `G` of a small elliptic curve `E`, together with the corresponding codomain
/// curve and point-evaluation formulas coming from Vélu's classical sums.
///
/// The core data are:
///
/// - a domain curve `E`
/// - a finite kernel subgroup `G`
/// - a codomain curve `E / G`
/// - the actual map `phi : E -> E / G`
pub struct VeluIsogeny<C: CurveModel> {
    pub(crate) domain: C,
    pub(crate) codomain: C,
    pub(crate) kernel: IsogenyKernel<C>,
    pub(crate) degree: usize,
}

impl<C: CurveModel + Clone> Clone for VeluIsogeny<C>
where
    C::Point: Clone,
{
    fn clone(&self) -> Self {
        Self {
            domain: self.domain.clone(),
            codomain: self.codomain.clone(),
            kernel: self.kernel.clone(),
            degree: self.degree,
        }
    }
}

impl<C: CurveModel> VeluIsogeny<C>
where
    C::Point: Clone + Eq + Hash,
{
    /// Builds a Vélu-isogeny scaffold from an explicit finite kernel.
    ///
    /// This is the natural entry point when the user already knows the whole
    /// subgroup `G` they want to collapse. The points are first validated as a
    /// genuine kernel subgroup of `domain`, and only then passed into the
    /// internal Vélu construction step.
    pub(crate) fn from_parts(domain: C, codomain: C, kernel: IsogenyKernel<C>) -> Self {
        let degree = kernel.degree();
        Self {
            domain,
            codomain,
            kernel,
            degree,
        }
    }

    /// Returns the validated finite kernel carried by the scaffold.
    ///
    /// This is the full subgroup `G`, including the identity point.
    pub fn kernel(&self) -> &IsogenyKernel<C> {
        &self.kernel
    }

    /// Returns the non-identity kernel points.
    ///
    /// Classical Vélu formulas usually sum over the non-zero points of the
    /// kernel. This helper relies on the [`crate::isogenies::IsogenyKernel`]
    /// storage invariant that the identity point is always stored first.
    pub fn kernel_nonzero_points(&self) -> &[C::Point] {
        self.kernel
            .points()
            .split_first()
            .map(|(_, tail)| tail)
            .unwrap_or(&[])
    }
}
