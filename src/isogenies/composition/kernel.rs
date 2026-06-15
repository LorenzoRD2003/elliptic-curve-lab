use crate::elliptic_curves::traits::{CurveModel, FiniteGroupCurveModel};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    composition::bridge::CompositionBridge, composition::evaluation::evaluate_composed_point,
    error::IsogenyError, traits::Isogeny,
};

pub(crate) fn compute_composed_kernel_points<Domain, Middle, Codomain, First, Second, Bridge>(
    first: &First,
    bridge: &Bridge,
    second: &Second,
) -> Result<Vec<Domain::Point>, IsogenyError>
where
    Domain: FiniteGroupCurveModel,
    Domain::BaseField: EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem>,
    Domain::Point: Clone + PartialEq,
    Middle: CurveModel,
    Codomain: CurveModel,
    Codomain::Point: PartialEq,
    First: Isogeny<Domain, Middle>,
    Second: Isogeny<Middle, Codomain>,
    Bridge: CompositionBridge<Middle>,
{
    let codomain_identity = second.codomain().identity();

    let kernel_points = first
        .domain()
        .points()
        .into_iter()
        .map(|point| {
            let image = evaluate_composed_point(first, bridge, second, &point)?;
            Ok((point, image == codomain_identity))
        })
        .collect::<Result<Vec<_>, IsogenyError>>()?
        .into_iter()
        .filter_map(|(point, maps_to_identity)| maps_to_identity.then_some(point))
        .collect::<Vec<_>>();

    Ok(kernel_points)
}
