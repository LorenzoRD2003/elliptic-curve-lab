use crate::elliptic_curves::traits::CurveModel;
use crate::isogenies::kernel::IsogenyKernel;

/// Reduced kernel data that can be exposed through explicit points.
#[derive(Clone, Debug)]
pub enum ReducedKernelDescription<C: CurveModel> {
    /// A validated explicit finite subgroup.
    RationalPointSubgroup(IsogenyKernel<C>),
    /// A finite reduced subgroup visible through explicit points, but not
    /// currently reified as an [`IsogenyKernel`].
    FiniteSubgroupSchemeVisibleAsPoints {
        /// Explicit points currently visible on the chosen base field.
        points: Vec<C::Point>,
        /// Degree contributed by the reduced subgroup.
        degree: usize,
    },
}

impl<C: CurveModel> ReducedKernelDescription<C> {
    /// Returns the reduced degree carried by this description.
    pub fn degree(&self) -> usize {
        match self {
            Self::RationalPointSubgroup(kernel) => kernel.degree(),
            Self::FiniteSubgroupSchemeVisibleAsPoints { degree, .. } => *degree,
        }
    }

    /// Returns the explicit points carried by this reduced description.
    pub fn points(&self) -> &[C::Point] {
        match self {
            Self::RationalPointSubgroup(kernel) => kernel.points(),
            Self::FiniteSubgroupSchemeVisibleAsPoints { points, .. } => points.as_slice(),
        }
    }
}
