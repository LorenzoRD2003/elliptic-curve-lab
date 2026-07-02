use crate::elliptic_curves::{
    CurveError,
    frobenius::{
        HasseInterval,
        hasse::search::{HasseMultipleSearchReport, HasseMultipleSearchStep},
    },
    traits::{BigScalarGroupCurveModel, HasseIntervalSearchCurveModel},
};
use num_bigint::BigUint;
use num_traits::ToPrimitive;

pub(crate) fn find_annihilating_multiple_in_interval_naive_report<
    C: HasseIntervalSearchCurveModel + ?Sized,
>(
    curve: &C,
    point: &C::Point,
    interval: HasseInterval,
) -> Result<HasseMultipleSearchReport<C::Point>, CurveError> {
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    let lower = interval.lower();
    let upper = interval.upper();
    let mut current = curve.mul_scalar_biguint(point, &lower)?;
    let mut steps = Vec::with_capacity(interval.candidate_count().to_usize().unwrap_or_default());
    let mut found = None;
    let mut candidate_multiple = lower.clone();

    while candidate_multiple <= upper {
        if candidate_multiple > lower {
            current = curve.add(&current, point)?;
        }

        steps.push(HasseMultipleSearchStep::new(
            candidate_multiple.clone(),
            current.clone(),
        ));

        if curve.is_identity(&current) {
            found = Some(candidate_multiple);
            break;
        }

        candidate_multiple += BigUint::from(1u8);
    }

    Ok(interval.search_report(found, steps))
}
