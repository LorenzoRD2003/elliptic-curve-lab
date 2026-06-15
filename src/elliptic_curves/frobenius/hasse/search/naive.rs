use num_bigint::BigUint;

use crate::elliptic_curves::{
    CurveError,
    frobenius::{
        HasseInterval,
        hasse::search::{HasseMultipleSearchReport, HasseMultipleSearchStep},
    },
    traits::{BigScalarGroupCurveModel, HasseIntervalSearchCurveModel},
};

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
    let mut current = curve.mul_scalar_biguint(point, &BigUint::from(lower))?;
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

    Ok(interval.search_report(found, steps))
}
