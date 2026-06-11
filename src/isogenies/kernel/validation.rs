use crate::elliptic_curves::traits::GroupCurveModel;
use crate::isogenies::{IsogenyError, IsogenyKernelError};
use std::collections::HashSet;
use std::hash::Hash;

pub(crate) fn validate_explicit_subgroup<C>(
    curve: &C,
    points: &HashSet<C::Point>,
) -> Result<(), IsogenyError>
where
    C: GroupCurveModel,
    C::Point: Clone + Eq + Hash,
{
    if points.is_empty() {
        return Err(IsogenyError::Kernel(IsogenyKernelError::EmptyKernel));
    }

    let identity = curve.identity();
    if !points.contains(&identity) {
        return Err(IsogenyError::Kernel(
            IsogenyKernelError::KernelDoesNotContainIdentity,
        ));
    }

    if points.iter().any(|point| !curve.contains(point)) {
        return Err(IsogenyError::Kernel(
            IsogenyKernelError::KernelPointNotOnCurve,
        ));
    }

    for point in points {
        let inverse = curve.neg(point);
        if !points.contains(&inverse) {
            return Err(IsogenyError::Kernel(
                IsogenyKernelError::KernelNotClosedUnderNegation,
            ));
        }
    }

    for left in points {
        for right in points {
            let sum = curve.add(left, right)?;
            if !points.contains(&sum) {
                return Err(IsogenyError::Kernel(
                    IsogenyKernelError::KernelNotClosedUnderAddition,
                ));
            }
        }
    }

    Ok(())
}
