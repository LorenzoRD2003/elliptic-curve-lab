use crate::elliptic_curves::isomorphisms::CurveIsomorphism;
use crate::elliptic_curves::traits::CurveModel;

use crate::isogenies::{Isogeny, IsogenyError, KernelDescription, ReducedKernelDescription};

/// Repackages an explicit curve isomorphism as a degree-one isogeny.
///
/// Mathematically, every curve isomorphism is in particular an isogeny of
/// degree `1`. This wrapper makes that perspective explicit inside the
/// isogeny module and supplies the small additional data that the current
/// educational [`Isogeny`] trait expects, namely the explicit rational kernel.
///
/// Since an isomorphism is injective, its rational kernel is just the identity
/// point.
pub struct IsomorphismIsogeny<Iso: CurveIsomorphism> {
    isomorphism: Iso,
    kernel_points: Vec<<Iso::Domain as CurveModel>::Point>,
}

impl<Iso: CurveIsomorphism> IsomorphismIsogeny<Iso> {
    /// Wraps an isomorphism as a degree-one isogeny.
    pub fn new(isomorphism: Iso) -> Self {
        let kernel_points = vec![isomorphism.domain().identity()];

        Self {
            isomorphism,
            kernel_points,
        }
    }

    /// Returns the wrapped isomorphism.
    pub fn isomorphism(&self) -> &Iso {
        &self.isomorphism
    }
}

impl<Iso: CurveIsomorphism> Isogeny<Iso::Domain, Iso::Codomain> for IsomorphismIsogeny<Iso> {
    fn domain(&self) -> &Iso::Domain {
        self.isomorphism.domain()
    }

    fn codomain(&self) -> &Iso::Codomain {
        self.isomorphism.codomain()
    }

    fn degree(&self) -> usize {
        1
    }

    fn evaluate(
        &self,
        point: &<Iso::Domain as CurveModel>::Point,
    ) -> Result<<Iso::Codomain as CurveModel>::Point, IsogenyError> {
        self.isomorphism.evaluate(point).map_err(Into::into)
    }

    fn kernel_description(&self) -> KernelDescription<Iso::Domain> {
        KernelDescription::Reduced(
            ReducedKernelDescription::FiniteSubgroupSchemeVisibleAsPoints {
                points: self.kernel_points.clone(),
                degree: 1,
            },
        )
    }
}
