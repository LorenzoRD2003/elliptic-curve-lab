use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, Zero};

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    short_weierstrass::rational_torsion::{
        RationalTorsionError, RationalTorsionGroup,
        enumeration::compare_affine_points,
        integral_model::{RationalIntegralModel, integral_rational_to_bigint},
        mazur::MAZUR_CYCLIC_ORDERS,
        reduction_mod_p::{
            good_prime::GoodReductionPrime,
            torsion_polynomial::{TorsionXPolynomial, TorsionXPolynomialError},
        },
    },
    traits::AffineCurveModel,
};
use crate::fields::Q;
use crate::numerics::exact_square_root;

/// Errors produced by the reduction/Hensel rational-point lift.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum ReducedTorsionLiftError {
    NoGoodReductionPrime,
    RationalTorsion(RationalTorsionError),
    TorsionPolynomial(TorsionXPolynomialError),
    Hensel(crate::numerics::hensel::HenselLiftError),
}

/// Verified rational torsion points recovered by the reduction/Hensel route.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ReducedTorsionLiftReport {
    #[cfg(test)]
    good_prime: GoodReductionPrime,
    group: RationalTorsionGroup,
    points: Vec<AffinePoint<Q>>,
}

impl ReducedTorsionLiftReport {
    /// Recovers rational torsion points on the integral model using the
    /// reduction/Hensel route.
    ///
    /// This stage stops at certified points on the integral companion model:
    /// it does not compare against Lutz-Nagell and does not transport points
    /// back to a non-integral source curve.
    ///
    /// Complexity: `Θ(M·p·n)`, where `M` is the fixed Mazur order-list length,
    /// `p` is the chosen good prime, and `n` is a coarse coefficient-count
    /// bound for the `x`-criteria.
    pub(super) fn from_integral_model(
        model: &RationalIntegralModel,
    ) -> Result<Self, ReducedTorsionLiftError> {
        let good_prime = GoodReductionPrime::first_for_integral_model(model)
            .ok_or(ReducedTorsionLiftError::NoGoodReductionPrime)?;
        debug_assert!(!good_prime.discriminant_mod_p().is_zero());
        let curve = model.curve();
        let mut verified = vec![(AffinePoint::infinity(), 1usize)];

        for order in MAZUR_CYCLIC_ORDERS {
            let polynomial = TorsionXPolynomial::from_integral_model(model, *order)?;
            let lift_report = polynomial.lift_roots_by_hensel(good_prime.prime())?;

            for x in lift_report.certified_roots() {
                for point in curve.points_for_x(&x)? {
                    if let Some(exact_order) = curve.exact_mazur_order(&point)? {
                        verified.push((point, exact_order));
                    }
                }
            }
        }

        let (points, point_orders): (Vec<_>, Vec<_>) =
            sort_and_dedup_verified_points(verified).into_iter().unzip();
        let group = RationalTorsionGroup::from_verified_point_orders(&point_orders)?;

        Ok(Self {
            #[cfg(test)]
            good_prime,
            group,
            points,
        })
    }

    #[cfg(test)]
    pub(super) fn good_prime(&self) -> GoodReductionPrime {
        self.good_prime
    }

    pub(super) fn group(&self) -> RationalTorsionGroup {
        self.group
    }

    pub(super) fn points(&self) -> &[AffinePoint<Q>] {
        &self.points
    }
}

impl ShortWeierstrassCurve<Q> {
    fn points_for_x(&self, x: &BigInt) -> Result<Vec<AffinePoint<Q>>, ReducedTorsionLiftError> {
        let b = integral_rational_to_bigint(self.b())?;
        let a = integral_rational_to_bigint(self.a())?;
        let rhs = x * x * x + a * x + b;
        if rhs.is_negative() {
            return Ok(Vec::new());
        }

        let Some(y_abs) = exact_square_root(rhs.magnitude()) else {
            return Ok(Vec::new());
        };
        let y_abs = BigInt::from(y_abs);
        let x = BigRational::from_integer(x.clone());

        if y_abs.is_zero() {
            return Ok(vec![self.point(x, BigRational::zero())?]);
        }

        let y = BigRational::from_integer(y_abs);
        Ok(vec![self.point(x.clone(), y.clone())?, self.point(x, -y)?])
    }
}

fn sort_and_dedup_verified_points(
    verified: Vec<(AffinePoint<Q>, usize)>,
) -> Vec<(AffinePoint<Q>, usize)> {
    let mut verified = verified;
    verified.sort_by(|(left, _), (right, _)| compare_affine_points(left, right));
    verified.dedup_by(|(left, _), (right, _)| left == right);
    verified
}

impl From<RationalTorsionError> for ReducedTorsionLiftError {
    fn from(error: RationalTorsionError) -> Self {
        Self::RationalTorsion(error)
    }
}

impl From<CurveError> for ReducedTorsionLiftError {
    fn from(error: CurveError) -> Self {
        Self::RationalTorsion(RationalTorsionError::from(error))
    }
}

impl From<TorsionXPolynomialError> for ReducedTorsionLiftError {
    fn from(error: TorsionXPolynomialError) -> Self {
        Self::TorsionPolynomial(error)
    }
}

impl From<crate::numerics::hensel::HenselLiftError> for ReducedTorsionLiftError {
    fn from(error: crate::numerics::hensel::HenselLiftError) -> Self {
        Self::Hensel(error)
    }
}

impl From<ReducedTorsionLiftError> for RationalTorsionError {
    fn from(error: ReducedTorsionLiftError) -> Self {
        let reason = match error {
            ReducedTorsionLiftError::NoGoodReductionPrime => "no good reduction prime was found",
            ReducedTorsionLiftError::RationalTorsion(error) => return error,
            ReducedTorsionLiftError::TorsionPolynomial(_) => {
                "a torsion x-criterion could not be constructed"
            }
            ReducedTorsionLiftError::Hensel(_) => "bounded simple Hensel lifting failed",
        };
        Self::ReductionHenselUnavailable { reason }
    }
}
