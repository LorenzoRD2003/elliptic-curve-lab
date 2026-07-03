use core::cmp::Ordering;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::elliptic_curves::{
    AffinePoint,
    short_weierstrass::rational_torsion::{
        RationalTorsionError, integral_model::RationalIntegralModel,
    },
    traits::AffineCurveModel,
};
use crate::fields::Q;
use crate::numerics::positive_divisors;
use crate::polynomials::IntegerPolynomial;

/// Bound used by the first rational-torsion route.
///
/// For `E(Q)_tors`, Mazur's theorem bounds every point order by the least
/// common multiple of the permitted cyclic orders:
/// `lcm(1, 2, ..., 10, 12) = 27720`.
pub(crate) const MAZUR_TORSION_EXPONENT_BOUND: usize = 27_720;

/// Cyclic orders permitted by Mazur's theorem over `Q`.
pub(crate) const MAZUR_CYCLIC_ORDERS: &[usize] = &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12];

/// Product-family parameters `m` for `ℤ/2ℤ × ℤ/2mℤ`.
pub(crate) const MAZUR_PRODUCT_PARAMETERS: &[usize] = &[1, 2, 3, 4];

/// Integral candidates produced by the Lutz-Nagell finite search.
///
/// For an integral model `E: y² = x³ + Ax + B`, Lutz-Nagell says every
/// non-identity torsion point in `E(Q)` has integral coordinates, and either
/// `y = 0` or `y² | Δ`. This report records the resulting finite candidate set
/// after exact curve-membership validation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct LutzNagellCandidateReport {
    discriminant: BigInt,
    points: Vec<AffinePoint<Q>>,
}

impl LutzNagellCandidateReport {
    fn new(discriminant: BigInt, points: Vec<AffinePoint<Q>>) -> Self {
        Self {
            discriminant,
            points,
        }
    }

    /// Builds the Lutz-Nagell candidate report from an integral model.
    ///
    /// Complexity: `Θ(factor(|Δ|) + Y·root_search + C log C)`, where `Y` is the
    /// number of possible `y`-values inspected and `C` is the number of
    /// validated candidates sorted for deterministic reporting.
    pub(crate) fn from_integral_model(
        model: &RationalIntegralModel,
    ) -> Result<Self, RationalTorsionError> {
        let curve = model.curve();
        let a = integral_rational_to_bigint(curve.a())?;
        let b = integral_rational_to_bigint(curve.b())?;
        let discriminant = integral_rational_to_bigint(&curve.discriminant())?;
        let mut candidates = vec![AffinePoint::infinity()];

        for y in lutz_nagell_y_values(&discriminant) {
            let y_squared = &y * &y;
            // solving x^3 + Ax + (B - y_i^2) = 0
            let polynomial = IntegerPolynomial::new(vec![
                &b - &y_squared,
                a.clone(),
                BigInt::zero(),
                BigInt::one(),
            ]);

            for x in polynomial.integer_roots_by_rational_root_test() {
                let point = curve.point(
                    BigRational::from_integer(x),
                    BigRational::from_integer(y.clone()),
                )?;
                candidates.push(point);
            }
        }

        sort_affine_points(&mut candidates);
        candidates.dedup();
        Ok(Self::new(discriminant, candidates))
    }

    pub(crate) fn discriminant(&self) -> &BigInt {
        &self.discriminant
    }

    pub(crate) fn points(&self) -> &[AffinePoint<Q>] {
        &self.points
    }

    pub(crate) fn candidate_count(&self) -> usize {
        self.points.len()
    }
}

/// Returns all integer `y`-values allowed by the Lutz-Nagell divisor test.
///
/// Complexity: `Θ(factor(|Δ|) + τ(|Δ|))`, where `τ(|Δ|)` is the number of
/// positive divisors inspected.
fn lutz_nagell_y_values(discriminant: &BigInt) -> Vec<BigInt> {
    let absolute_discriminant = discriminant.magnitude();
    let mut values = vec![BigInt::zero()];

    for divisor in positive_divisors(absolute_discriminant) {
        if (&divisor * &divisor).is_zero() {
            continue;
        }
        if (absolute_discriminant % (&divisor * &divisor)).is_zero() {
            let positive = BigInt::from(divisor);
            values.push(-&positive);
            values.push(positive);
        }
    }

    values.sort();
    values.dedup();
    values
}

fn integral_rational_to_bigint(value: &BigRational) -> Result<BigInt, RationalTorsionError> {
    if !value.is_integer() {
        return Err(RationalTorsionError::IntegralModelUnavailable);
    }

    Ok(value.to_integer())
}

fn sort_affine_points(points: &mut [AffinePoint<Q>]) {
    points.sort_by(compare_affine_points);
}

fn compare_affine_points(left: &AffinePoint<Q>, right: &AffinePoint<Q>) -> Ordering {
    match (left, right) {
        (AffinePoint::Infinity, AffinePoint::Infinity) => Ordering::Equal,
        (AffinePoint::Infinity, AffinePoint::Finite { .. }) => Ordering::Less,
        (AffinePoint::Finite { .. }, AffinePoint::Infinity) => Ordering::Greater,
        (
            AffinePoint::Finite {
                x: left_x,
                y: left_y,
            },
            AffinePoint::Finite {
                x: right_x,
                y: right_y,
            },
        ) => left_x.cmp(right_x).then_with(|| left_y.cmp(right_y)),
    }
}
