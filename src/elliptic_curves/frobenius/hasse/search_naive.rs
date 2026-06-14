use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::frobenius::{
    HasseInterval, HasseMultipleSearchReport, HasseMultipleSearchStep, hasse_multiple_search_report,
};
use crate::elliptic_curves::order_from_multiple::mul_scalar_biguint;
use crate::elliptic_curves::traits::GroupCurveModel;
use num_bigint::BigUint;

/// Searches one already-chosen interval from left to right until `[M]P = O`
/// is found or the interval is exhausted.
///
/// Complexity: one `BigUint` scalar multiplication to build `[L]P`, then
/// `Θ(|H|)` group additions, where `|H|` is the number of integer candidates
/// in the supplied interval.
pub(crate) fn find_annihilating_multiple_in_interval_naive_report<C>(
    curve: &C,
    point: &C::Point,
    interval: HasseInterval,
) -> Result<HasseMultipleSearchReport<C::Point>, CurveError>
where
    C: GroupCurveModel + ?Sized,
    C::Point: Clone,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    let lower = interval.lower();
    let upper = interval.upper();
    let mut current = mul_scalar_biguint(curve, point, &BigUint::from(lower))?;
    let mut steps = Vec::with_capacity(interval.candidate_count() as usize);
    let mut found = None;

    for candidate_multiple in lower..=upper {
        if candidate_multiple > lower {
            current = curve.add(&current, point)?;
        }

        steps.push(HasseMultipleSearchStep::new(
            candidate_multiple,
            current.clone(),
        ));

        if curve.is_identity(&current) {
            found = Some(candidate_multiple);
            break;
        }
    }

    Ok(hasse_multiple_search_report(
        interval.q(),
        interval,
        found,
        steps,
    ))
}
