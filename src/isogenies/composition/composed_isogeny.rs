use core::marker::PhantomData;

use crate::elliptic_curves::traits::{CurveIsomorphism, CurveModel, FiniteGroupCurveModel};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    composition::{
        bridge::CompositionBridge, evaluation::evaluate_composed_point,
        kernel::compute_composed_kernel_points,
    },
    error::IsogenyError,
    kernel::{KernelDescription, ReducedKernelDescription},
    traits::Isogeny,
};

/// Formal composition `second ∘ first` of two explicit isogenies.
///
/// This educational scaffold models the usual right-to-left composition:
///
/// - `first : E -> E'`
/// - `second : E' -> E''`
/// - `second ∘ first : E -> E''`
///
/// The field names intentionally use `first` and `second` instead of
/// `left` and `right`, since composition is read right-to-left and those names
/// are often easier to follow in educational code and prose.
pub struct ComposedIsogeny<I, J, Domain: CurveModel, Middle, Codomain, Bridge = ()> {
    first: I,
    bridge: Bridge,
    second: J,
    kernel_points: Vec<Domain::Point>,
    marker: PhantomData<(Domain, Middle, Codomain)>,
}

impl<I, J, Domain, Middle, Codomain> ComposedIsogeny<I, J, Domain, Middle, Codomain, ()>
where
    Domain: FiniteGroupCurveModel,
    Domain::BaseField: EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem>,
    Domain::Point: Clone + PartialEq,
    Middle: CurveModel + PartialEq,
    Codomain: CurveModel,
    Codomain::Point: PartialEq,
    I: Isogeny<Domain, Middle>,
    J: Isogeny<Middle, Codomain>,
{
    /// Returns the first map in the composition.
    pub fn first(&self) -> &I {
        &self.first
    }

    /// Returns the second map in the composition.
    pub fn second(&self) -> &J {
        &self.second
    }

    /// Builds the formal composition `second ∘ first` under strict middle-curve equality.
    pub fn new_strict(first: I, second: J) -> Result<Self, IsogenyError> {
        ().validate_bridge(first.codomain(), second.domain())?;
        let bridge = ();
        let kernel_points = compute_composed_kernel_points(&first, &bridge, &second)?;

        Ok(Self {
            first,
            bridge,
            second,
            kernel_points,
            marker: PhantomData,
        })
    }
}

impl<I, J, Domain, Middle, Codomain, Bridge> ComposedIsogeny<I, J, Domain, Middle, Codomain, Bridge>
where
    Domain: FiniteGroupCurveModel,
    Domain::BaseField: EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem>,
    Domain::Point: Clone + PartialEq,
    Middle: CurveModel + Clone + PartialEq,
    Codomain: CurveModel,
    Codomain::Point: PartialEq,
    I: Isogeny<Domain, Middle>,
    J: Isogeny<Middle, Codomain>,
    Bridge: CurveIsomorphism<Domain = Middle, Codomain = Middle>,
{
    /// Builds the bridged composition `second ∘ bridge ∘ first`.
    pub fn new_up_to_isomorphism(
        first: I,
        bridge: Bridge,
        second: J,
    ) -> Result<Self, IsogenyError> {
        bridge.validate_bridge(first.codomain(), second.domain())?;
        let kernel_points = compute_composed_kernel_points(&first, &bridge, &second)?;

        Ok(Self {
            first,
            bridge,
            second,
            kernel_points,
            marker: PhantomData,
        })
    }
}

impl<I, J, Domain, Middle, Codomain, Bridge> Isogeny<Domain, Codomain>
    for ComposedIsogeny<I, J, Domain, Middle, Codomain, Bridge>
where
    Domain: FiniteGroupCurveModel,
    Domain::BaseField: EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem>,
    Domain::Point: Clone + PartialEq,
    Middle: CurveModel + Clone + PartialEq,
    Codomain: CurveModel,
    Codomain::Point: PartialEq,
    I: Isogeny<Domain, Middle>,
    J: Isogeny<Middle, Codomain>,
    Bridge: CompositionBridge<Middle>,
{
    fn domain(&self) -> &Domain {
        self.first.domain()
    }

    fn codomain(&self) -> &Codomain {
        self.second.codomain()
    }

    fn degree(&self) -> usize {
        self.first.degree() * self.second.degree()
    }

    fn evaluate(&self, point: &Domain::Point) -> Result<Codomain::Point, IsogenyError> {
        evaluate_composed_point(&self.first, &self.bridge, &self.second, point)
    }

    fn kernel_description(&self) -> KernelDescription<Domain> {
        KernelDescription::Reduced(
            ReducedKernelDescription::FiniteSubgroupSchemeVisibleAsPoints {
                points: self.kernel_points.clone(),
                degree: self.kernel_points.len(),
            },
        )
    }
}
