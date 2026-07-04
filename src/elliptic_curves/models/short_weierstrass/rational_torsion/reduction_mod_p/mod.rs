//! Internal reduction-mod-`p` helpers for the rational-torsion route.
//!
//! This module deliberately does not introduce a `Field` implementation.  The
//! stage-5 algorithm needs a small runtime prime chosen from the integral
//! model, while `fields::Fp` is intentionally static and type-level.

#[cfg(test)]
mod formula_ops;
mod good_prime;
mod hensel_lift;
mod rational_points;
#[cfg(test)]
mod reduced_curve;
mod small_prime_field;
mod torsion_polynomial;

#[cfg(test)]
mod tests;

use crate::elliptic_curves::{
    AffinePoint,
    short_weierstrass::rational_torsion::{
        RationalTorsionError, RationalTorsionGroup, integral_model::RationalIntegralModel,
    },
};
use crate::fields::Q;

pub(super) fn rational_points_from_integral_model(
    model: &RationalIntegralModel,
) -> Result<(RationalTorsionGroup, Vec<AffinePoint<Q>>), RationalTorsionError> {
    let report = rational_points::ReducedTorsionLiftReport::from_integral_model(model)?;
    Ok((report.group(), report.points().to_vec()))
}
