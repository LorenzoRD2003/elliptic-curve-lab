use crate::fields::{EnumerableFiniteField, SqrtField};

use super::{CurveError, FiniteGroupCurveModel};

fn proper_divisors_from_prime_factors(n: usize) -> Vec<usize> {
    let mut divisors = Vec::new();
    let mut remaining = n;
    let mut factor = 2;

    while factor * factor <= remaining {
        if remaining.is_multiple_of(factor) {
            divisors.push(n / factor);
            while remaining.is_multiple_of(factor) {
                remaining /= factor;
            }
        }

        factor += if factor == 2 { 1 } else { 2 };
    }

    if remaining > 1 && n > 1 {
        divisors.push(n / remaining);
    }

    divisors
}

/// Returns whether `point` has exact order `n`.
///
/// This is the generic small-group torsion helper that the codebase uses when
/// the full curve group can be enumerated honestly.
///
/// The logic follows the classical exact-order criterion:
///
/// - first check `[n]P = O`
/// - then rule out every quotient `n / p` where `p` is a prime divisor of `n`
///
/// This is enough because a point has exact order `n` exactly when it is
/// killed by `[n]` but by no `n / p` for prime `p | n`.
pub fn point_has_exact_order<C>(curve: &C, point: &C::Point, n: usize) -> Result<bool, CurveError>
where
    C: FiniteGroupCurveModel,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
{
    if n == 0 {
        return Err(CurveError::InvalidTorsionOrder { order: n });
    }

    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    if curve.is_identity(point) {
        return Ok(n == 1);
    }

    let n_multiple = curve.mul_scalar(point, n as u64)?;
    if !curve.is_identity(&n_multiple) {
        return Ok(false);
    }

    for divisor in proper_divisors_from_prime_factors(n) {
        let divisor_multiple = curve.mul_scalar(point, divisor as u64)?;
        if curve.is_identity(&divisor_multiple) {
            return Ok(false);
        }
    }

    Ok(true)
}

/// Returns the non-identity rational points of exact order `order`.
///
/// This helper reuses [`point_has_exact_order`] over the fully enumerated
/// rational point set, so it is honest for both prime and composite orders in
/// the current small finite educational setting.
pub fn points_of_exact_order<C>(curve: &C, order: usize) -> Result<Vec<C::Point>, CurveError>
where
    C: FiniteGroupCurveModel,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
{
    if order == 0 {
        return Err(CurveError::InvalidTorsionOrder { order });
    }

    let mut points = Vec::new();
    for point in curve.points() {
        if curve.is_identity(&point) {
            continue;
        }

        if point_has_exact_order(curve, &point, order)? {
            points.push(point);
        }
    }

    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::{point_has_exact_order, points_of_exact_order, proper_divisors_from_prime_factors};
    use crate::elliptic_curves::{
        AffineCurveModel, AffinePoint, CurveError, CurveModel, ShortWeierstrassCurve,
    };
    use crate::fields::{Field, Fp};

    type F5 = Fp<5>;
    type F7 = Fp<7>;
    type F41 = Fp<41>;
    type Curve5 = ShortWeierstrassCurve<F5>;
    type Curve7 = ShortWeierstrassCurve<F7>;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn f5_noncyclic_curve() -> Curve5 {
        Curve5::new(F5::from_i64(-1), F5::zero()).expect("valid curve")
    }

    fn f7_curve() -> Curve7 {
        Curve7::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
    }

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn proper_divisors_from_prime_factors_uses_only_prime_factor_checks() {
        assert_eq!(proper_divisors_from_prime_factors(2), vec![1]);
        assert_eq!(proper_divisors_from_prime_factors(12), vec![6, 4]);
        assert_eq!(proper_divisors_from_prime_factors(27), vec![9]);
    }

    #[test]
    fn point_has_exact_order_rejects_zero_degree() {
        let curve = f41_curve();
        let point = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");

        assert_eq!(
            point_has_exact_order(&curve, &point, 0),
            Err(CurveError::InvalidTorsionOrder { order: 0 })
        );
    }

    #[test]
    fn point_has_exact_order_rejects_points_outside_the_curve() {
        let curve = f41_curve();
        let off_curve = AffinePoint::<F41>::new(F41::from_i64(2), F41::from_i64(2));

        assert_eq!(
            point_has_exact_order(&curve, &off_curve, 2),
            Err(CurveError::PointNotOnCurve)
        );
    }

    #[test]
    fn point_has_exact_order_accepts_exact_order_two_point() {
        let curve = f41_curve();
        let point = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");

        assert_eq!(point_has_exact_order(&curve, &point, 2), Ok(true));
    }

    #[test]
    fn point_has_exact_order_rejects_points_killed_by_a_prime_factor_quotient() {
        let curve = f41_curve();
        let point = curve
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("sample point should lie on the curve");

        assert_eq!(point_has_exact_order(&curve, &point, 4), Ok(false));
    }

    #[test]
    fn points_of_exact_order_rejects_zero_degree() {
        let curve = f41_curve();

        assert_eq!(
            points_of_exact_order(&curve, 0),
            Err(CurveError::InvalidTorsionOrder { order: 0 })
        );
    }

    #[test]
    fn points_of_exact_order_returns_non_identity_points_of_exact_order() {
        let curve = f41_curve();
        let points = points_of_exact_order(&curve, 2).expect("degree two should be valid");

        assert_eq!(
            points,
            vec![curve.point(F41::from_i64(40), F41::from_i64(0)).unwrap()]
        );
        assert!(points.iter().all(|point| !curve.is_identity(point)));
    }

    #[test]
    fn finds_nonzero_two_torsion_points() {
        let curve = f5_noncyclic_curve();
        let points = points_of_exact_order(&curve, 2).expect("degree two should be valid");

        assert_eq!(
            points,
            vec![
                curve.point(F5::from_i64(0), F5::from_i64(0)).unwrap(),
                curve.point(F5::from_i64(1), F5::from_i64(0)).unwrap(),
                curve.point(F5::from_i64(4), F5::from_i64(0)).unwrap(),
            ]
        );
    }

    #[test]
    fn identity_is_not_reported_as_nontrivial_exact_order_point() {
        let curve = f5_noncyclic_curve();
        let identity = curve.identity();

        assert_eq!(point_has_exact_order(&curve, &identity, 2), Ok(false));

        let points = points_of_exact_order(&curve, 2).expect("degree two should be valid");
        assert!(!points.contains(&identity));
    }

    #[test]
    fn points_of_exact_order_matches_known_small_order_three_example() {
        let curve = f7_curve();
        let points = points_of_exact_order(&curve, 3).expect("degree three should be valid");

        assert_eq!(
            points,
            vec![
                curve.point(F7::from_i64(3), F7::from_i64(1)).unwrap(),
                curve.point(F7::from_i64(3), F7::from_i64(6)).unwrap(),
            ]
        );
    }
}
