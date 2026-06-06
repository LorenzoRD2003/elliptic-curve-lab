use std::collections::HashSet;
use std::hash::Hash;

use crate::elliptic_curves::traits::{CurveModel, FiniteGroupCurveModel};
use crate::fields::{EnumerableFiniteField, SqrtField};

use crate::isogenies::{Isogeny, IsogenyError, IsogenyKernel};

use crate::isogenies::velu::SupportsVeluConstruction;

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
    pub(super) domain: C,
    pub(super) codomain: C,
    pub(super) kernel: IsogenyKernel<C>,
    pub(super) degree: usize,
}

impl<C> Clone for VeluIsogeny<C>
where
    C: CurveModel + Clone,
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

#[allow(private_bounds)]
impl<C> VeluIsogeny<C>
where
    C: SupportsVeluConstruction,
    C::Point: Clone + Eq + Hash,
{
    /// Builds a Vélu-isogeny scaffold from an explicit finite kernel.
    ///
    /// This is the natural entry point when the user already knows the whole
    /// subgroup `G` they want to collapse. The points are first validated as a
    /// genuine kernel subgroup of `domain`, and only then passed into the
    /// internal Vélu construction step.
    pub fn from_points(domain: C, kernel_points: HashSet<C::Point>) -> Result<Self, IsogenyError> {
        let kernel = IsogenyKernel::new(&domain, kernel_points)?;
        Self::from_kernel(domain, kernel)
    }

    /// Builds a Vélu-isogeny scaffold from a cyclic kernel generator.
    ///
    /// This is the most convenient entry point for first-pass experiments:
    /// start from a torsion point `P`, build the cyclic subgroup `<P>`, and
    /// then hand that validated subgroup to the internal Vélu construction.
    pub fn from_generator(domain: C, generator: C::Point) -> Result<Self, IsogenyError>
    where
        C: FiniteGroupCurveModel,
        C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    {
        let kernel = IsogenyKernel::cyclic(&domain, &generator)?;
        Self::from_kernel(domain, kernel)
    }

    /// This method starts from a curve `E` and a finite subgroup `G`, then
    /// computes the codomain curve and explicit map data directly from Vélu's formulas.
    pub(super) fn from_kernel(domain: C, kernel: IsogenyKernel<C>) -> Result<Self, IsogenyError> {
        let codomain = C::velu_codomain_from_kernel(&domain, &kernel)?;
        let degree = kernel.degree();

        Ok(Self {
            domain,
            codomain,
            kernel,
            degree,
        })
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
    /// kernel. This helper relies on the [`IsogenyKernel`] storage invariant
    /// that the identity point is always stored first, so the non-zero kernel
    /// points are exactly the tail of the kernel slice.
    pub fn kernel_nonzero_points(&self) -> &[C::Point] {
        self.kernel
            .points()
            .split_first()
            .map(|(_, tail)| tail)
            .unwrap_or(&[])
    }
}

#[allow(private_bounds)]
impl<C> Isogeny<C, C> for VeluIsogeny<C>
where
    C: SupportsVeluConstruction,
    C::Point: Clone + Eq + Hash,
{
    fn domain(&self) -> &C {
        &self.domain
    }

    fn codomain(&self) -> &C {
        &self.codomain
    }

    fn degree(&self) -> usize {
        self.degree
    }

    fn evaluate(&self, point: &C::Point) -> Result<C::Point, IsogenyError> {
        if !self.domain.contains(point) {
            return Err(IsogenyError::Curve(
                crate::elliptic_curves::CurveError::PointNotOnCurve,
            ));
        }

        if self.kernel.contains(point) {
            return Ok(self.codomain.identity());
        }

        let image = C::velu_evaluate_non_kernel_point(self, point)?;
        debug_assert!(self.codomain.contains(&image));
        Ok(image)
    }

    fn kernel_points(&self) -> &[C::Point] {
        self.kernel.points()
    }
}
