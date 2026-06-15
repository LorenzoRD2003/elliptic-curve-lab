use crate::elliptic_curves::{CurveError, traits::CurveModel};
use crate::isogenies::{
    composition::bridge::CompositionBridge,
    error::{IsogenyError, IsogenyVerificationError},
    traits::Isogeny,
};

pub(crate) fn evaluate_composed_point<Domain, Middle, Codomain, First, Second, Bridge>(
    first: &First,
    bridge: &Bridge,
    second: &Second,
    point: &Domain::Point,
) -> Result<Codomain::Point, IsogenyError>
where
    Domain: CurveModel,
    Middle: CurveModel,
    Codomain: CurveModel,
    First: Isogeny<Domain, Middle>,
    Second: Isogeny<Middle, Codomain>,
    Bridge: CompositionBridge<Middle>,
{
    if !first.domain().contains(point) {
        return Err(CurveError::PointNotOnCurve.into());
    }

    let middle = first.evaluate(point)?;
    if !first.codomain().contains(&middle) {
        return Err(IsogenyError::Verification(
            IsogenyVerificationError::ImagePointNotOnCodomain,
        ));
    }

    let bridged = bridge.transport(middle)?;
    if !second.domain().contains(&bridged) {
        return Err(IsogenyError::Verification(
            IsogenyVerificationError::ImagePointNotOnCodomain,
        ));
    }

    let image = second.evaluate(&bridged)?;
    if !second.codomain().contains(&image) {
        return Err(IsogenyError::Verification(
            IsogenyVerificationError::ImagePointNotOnCodomain,
        ));
    }

    Ok(image)
}
