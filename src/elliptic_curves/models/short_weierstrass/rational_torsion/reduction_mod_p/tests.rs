use num_bigint::{BigInt, BigUint};
use num_rational::BigRational;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::rational_torsion::{
        RationalTorsionGroupShape,
        integral_model::RationalIntegralModel,
        mazur::MAZUR_CYCLIC_ORDERS,
        reduction_mod_p::{
            good_prime::GoodReductionPrime,
            hensel_lift::{TorsionXSeedLiftOutcome, lift_polynomial_roots_by_hensel},
            rational_points::ReducedTorsionLiftReport,
            reduced_curve::{ReducedPoint, ReducedShortWeierstrassCurve},
            small_prime_field::{ReductionPrime, ReductionPrimeError},
            torsion_polynomial::{
                TorsionXPolynomial, TorsionXPolynomialError, TorsionXPolynomialSource,
            },
        },
    },
    traits::CurveModel,
};
use crate::fields::Q;
use crate::polynomials::IntegerPolynomial;

fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn integral_model(a: i64, b: i64) -> RationalIntegralModel {
    let curve = ShortWeierstrassCurve::<Q>::new(q(a, 1), q(b, 1))
        .expect("test curve should be nonsingular");
    RationalIntegralModel::from_curve(curve).expect("test curve should already be integral")
}

fn rational_point(x: i64, y: i64) -> AffinePoint<Q> {
    AffinePoint::new(q(x, 1), q(y, 1))
}

#[test]
fn reduction_prime_constructor_accepts_only_primes() {
    assert_eq!(ReductionPrime::new(0), Err(ReductionPrimeError::Zero));
    assert_eq!(ReductionPrime::new(1), Err(ReductionPrimeError::One));
    assert_eq!(
        ReductionPrime::new(21),
        Err(ReductionPrimeError::Composite { modulus: 21 })
    );
    assert_eq!(ReductionPrime::new(17).expect("17 is prime").modulus(), 17);
}

#[test]
fn reduction_prime_reduces_bigints_to_canonical_residues() {
    let prime = ReductionPrime::new(17).expect("17 is prime");

    assert_eq!(prime.reduce_bigint(&BigInt::from(0)).representative(), 0);
    assert_eq!(prime.reduce_bigint(&BigInt::from(35)).representative(), 1);
    assert_eq!(prime.reduce_bigint(&BigInt::from(-1)).representative(), 16);
    assert_eq!(prime.reduce_bigint(&BigInt::from(-34)).representative(), 0);
}

#[test]
fn reduction_prime_arithmetic_keeps_canonical_representatives() {
    let prime = ReductionPrime::new(17).expect("17 is prime");
    let five = prime.reduce_u64(5);
    let thirteen = prime.reduce_u64(13);

    assert_eq!(prime.add(five, thirteen).representative(), 1);
    assert_eq!(prime.sub(five, thirteen).representative(), 9);
    assert_eq!(prime.mul(five, thirteen).representative(), 14);
    assert_eq!(prime.neg(five).representative(), 12);
    assert_eq!(prime.pow(five, 3).representative(), 6);
}

#[test]
fn reduction_prime_inverts_nonzero_residues() {
    let prime = ReductionPrime::new(17).expect("17 is prime");
    let five = prime.reduce_u64(5);
    let inverse = prime.inv(five).expect("5 is non-zero modulo 17");

    assert_eq!(inverse.representative(), 7);
    assert_eq!(prime.mul(five, inverse), prime.one());
    assert_eq!(prime.inv(prime.zero()), None);
}

#[test]
fn reduction_prime_enumerates_every_residue_once() {
    let prime = ReductionPrime::new(11).expect("11 is prime");
    let residues = prime
        .residues()
        .map(|residue| residue.representative())
        .collect::<Vec<_>>();

    assert_eq!(residues, (0..11).collect::<Vec<_>>());
}

#[test]
fn good_reduction_prime_starts_at_eleven() {
    let model = integral_model(-1, 0);
    let good_prime = GoodReductionPrime::first_for_integral_model(&model)
        .expect("y² = x³ - x has good reduction at 11");

    assert_eq!(good_prime.prime().modulus(), 11);
}

#[test]
fn good_reduction_prime_skips_primes_dividing_discriminant() {
    let model = integral_model(2, 3);
    let good_prime = GoodReductionPrime::first_for_integral_model(&model)
        .expect("a later small prime should give good reduction");

    assert_eq!(good_prime.prime().modulus(), 13);
}

#[test]
fn good_reduction_prime_records_nonzero_discriminant_residue() {
    let model = integral_model(1, 1);
    let good_prime = GoodReductionPrime::first_for_integral_model(&model)
        .expect("the fixture should have a good small prime");

    assert_eq!(good_prime.prime().modulus(), 11);
    assert_ne!(good_prime.discriminant_mod_p(), good_prime.prime().zero());
}

#[test]
fn reduced_curve_stores_coefficients_modulo_the_good_prime() {
    let model = integral_model(-1, 0);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");

    assert_eq!(reduced.prime().modulus(), 11);
    assert_eq!(reduced.a().representative(), 10);
    assert_eq!(reduced.b().representative(), 0);
    assert_ne!(
        reduced.good_prime().discriminant_mod_p(),
        reduced.prime().zero()
    );
}

#[test]
fn reduced_curve_enumerates_points_on_y_squared_equals_x_cubed_minus_x() {
    let model = integral_model(-1, 0);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");
    let points = reduced.points();

    assert_eq!(points.len(), 12);
    assert_eq!(points[0], ReducedPoint::Infinity);
    assert!(points.iter().copied().all(|point| reduced.contains(point)));
    assert!(points.contains(&ReducedPoint::Finite {
        x: reduced.prime().reduce_u64(4),
        y: reduced.prime().reduce_u64(7),
    }));
}

#[test]
fn reduced_curve_uses_the_later_good_prime_when_eleven_is_bad() {
    let model = integral_model(2, 3);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");

    assert_eq!(reduced.prime().modulus(), 13);
    assert_eq!(reduced.a().representative(), 2);
    assert_eq!(reduced.b().representative(), 3);
    assert_eq!(reduced.points().len(), 18);
}

#[test]
fn reduced_curve_group_law_detects_small_orders() {
    let model = integral_model(-1, 0);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");
    let prime = reduced.prime();
    let point = ReducedPoint::Finite {
        x: prime.reduce_u64(4),
        y: prime.reduce_u64(4),
    };
    let negated = ReducedPoint::Finite {
        x: prime.reduce_u64(4),
        y: prime.reduce_u64(7),
    };

    assert_eq!(reduced.neg(point), negated);
    assert_eq!(reduced.add(point, negated), ReducedPoint::Infinity);
    assert_eq!(reduced.mul_scalar(point, 2), negated);
    assert_eq!(reduced.mul_scalar(point, 3), ReducedPoint::Infinity);
    assert_eq!(reduced.mazur_order(point), Some(3));
}

#[test]
fn reduced_curve_filters_points_by_mazur_permitted_orders() {
    let model = integral_model(2, 3);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");
    let filtered = reduced.mazur_order_points();

    assert_eq!(reduced.points().len(), 18);
    assert_eq!(filtered.len(), 12);
    assert!(filtered.iter().all(|point| reduced.contains(point.point())));
    assert!(filtered.iter().all(|point| point.order() <= 12));
    assert_eq!(
        reduced.mazur_order(ReducedPoint::Finite {
            x: reduced.prime().reduce_u64(4),
            y: reduced.prime().reduce_u64(6),
        }),
        None
    );
}

#[test]
fn torsion_x_polynomial_for_order_two_is_the_defining_cubic() {
    let model = integral_model(-1, 0);
    let polynomial = TorsionXPolynomial::from_integral_model(&model, 2)
        .expect("order-two criterion should be available");

    assert_eq!(polynomial.order(), 2);
    assert_eq!(
        polynomial.source(),
        TorsionXPolynomialSource::TwoTorsionCubic
    );
    assert_eq!(
        polynomial.polynomial().to_dense_coefficients(),
        vec![
            BigInt::from(0),
            BigInt::from(-1),
            BigInt::from(0),
            BigInt::from(1),
        ]
    );
}

#[test]
fn torsion_x_polynomial_reuses_odd_division_polynomial() {
    let model = integral_model(-1, 0);
    let polynomial = TorsionXPolynomial::from_integral_model(&model, 3)
        .expect("order-three criterion should be available");

    assert_eq!(polynomial.order(), 3);
    assert_eq!(
        polynomial.source(),
        TorsionXPolynomialSource::OddDivisionPolynomial
    );
    assert_eq!(
        polynomial.polynomial().to_dense_coefficients(),
        vec![
            BigInt::from(-1),
            BigInt::from(0),
            BigInt::from(-6),
            BigInt::from(0),
            BigInt::from(3),
        ]
    );
}

#[test]
fn torsion_x_polynomial_reuses_even_division_polynomial_factor() {
    let model = integral_model(-1, 0);
    let polynomial = TorsionXPolynomial::from_integral_model(&model, 4)
        .expect("order-four criterion should be available");

    assert_eq!(polynomial.order(), 4);
    assert_eq!(
        polynomial.source(),
        TorsionXPolynomialSource::EvenDivisionPolynomialOverPsi2
    );
    assert_eq!(
        polynomial.polynomial().to_dense_coefficients(),
        vec![
            BigInt::from(1),
            BigInt::from(0),
            BigInt::from(-5),
            BigInt::from(0),
            BigInt::from(-5),
            BigInt::from(0),
            BigInt::from(1),
        ]
    );
}

#[test]
fn torsion_x_polynomial_exists_for_every_mazur_cyclic_order() {
    let model = integral_model(1, 1);

    for order in MAZUR_CYCLIC_ORDERS {
        let polynomial = TorsionXPolynomial::from_integral_model(&model, *order)
            .expect("every Mazur cyclic order should have an x-criterion");

        assert_eq!(polynomial.order(), *order);
        assert!(!polynomial.polynomial().is_zero());
    }
}

#[test]
fn torsion_x_polynomial_rejects_orders_outside_mazur_list() {
    let model = integral_model(1, 1);

    assert_eq!(
        TorsionXPolynomial::from_integral_model(&model, 1),
        Err(TorsionXPolynomialError::UnsupportedOrder { order: 1 })
    );
    assert_eq!(
        TorsionXPolynomial::from_integral_model(&model, 11),
        Err(TorsionXPolynomialError::UnsupportedOrder { order: 11 })
    );
}

#[test]
fn torsion_x_polynomial_finds_order_two_roots_modulo_good_prime() {
    let model = integral_model(-1, 0);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");
    let polynomial = TorsionXPolynomial::from_integral_model(&model, 2)
        .expect("order-two criterion should be available");
    let roots = polynomial
        .roots_mod_prime(reduced.prime())
        .into_iter()
        .map(|root| root.representative())
        .collect::<Vec<_>>();

    assert_eq!(roots, vec![0, 1, 10]);
}

#[test]
fn torsion_x_polynomial_vanishes_on_reduced_mazur_order_points() {
    let model = integral_model(2, 3);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");

    for point in reduced.mazur_order_points() {
        let ReducedPoint::Finite { x, .. } = point.point() else {
            continue;
        };
        let polynomial = TorsionXPolynomial::from_integral_model(&model, point.order())
            .expect("Mazur order criterion should be available");

        assert_eq!(
            polynomial.evaluate_mod_prime(reduced.prime(), x),
            reduced.prime().zero()
        );
    }
}

#[test]
fn torsion_x_hensel_lift_certifies_order_two_integer_roots() {
    let model = integral_model(-1, 0);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");
    let polynomial = TorsionXPolynomial::from_integral_model(&model, 2)
        .expect("order-two criterion should be available");
    let report = polynomial
        .lift_roots_by_hensel(reduced.prime())
        .expect("simple order-two seeds should lift cleanly");

    assert_eq!(report.order(), 2);
    assert_eq!(report.source(), TorsionXPolynomialSource::TwoTorsionCubic);
    assert_eq!(report.prime().modulus(), 11);
    assert_eq!(report.root_bound(), &BigUint::from(2u8));
    assert_eq!(report.seeds().len(), 3);
    assert_eq!(
        report.certified_roots(),
        vec![BigInt::from(-1), BigInt::from(0), BigInt::from(1)]
    );
    assert_eq!(report.singular_seed_count(), 0);
    assert_eq!(report.uncertified_seed_count(), 0);
}

#[test]
fn torsion_x_hensel_lift_runs_even_psi_m_over_psi_two_route() {
    let model = integral_model(0, 1);
    let reduced = ReducedShortWeierstrassCurve::from_integral_model(&model)
        .expect("fixture should have a good reduction prime");
    let polynomial = TorsionXPolynomial::from_integral_model(&model, 4)
        .expect("order-four criterion should be available");
    let report = polynomial
        .lift_roots_by_hensel(reduced.prime())
        .expect("order-four seeds should be reportable");

    assert_eq!(report.order(), 4);
    assert_eq!(
        report.source(),
        TorsionXPolynomialSource::EvenDivisionPolynomialOverPsi2
    );
    assert_eq!(report.prime().modulus(), 11);
    assert!(!report.seeds().is_empty());
}

#[test]
fn torsion_x_hensel_lift_reports_singular_seeds() {
    let prime = ReductionPrime::new(5).expect("5 is prime");
    let polynomial =
        IntegerPolynomial::new(vec![BigInt::from(0), BigInt::from(0), BigInt::from(1)]);
    let report = lift_polynomial_roots_by_hensel(
        3,
        TorsionXPolynomialSource::OddDivisionPolynomial,
        &polynomial,
        prime,
        vec![prime.zero()],
    )
    .expect("singular seeds should be non-fatal");

    assert_eq!(report.singular_seed_count(), 1);
    assert_eq!(report.uncertified_seed_count(), 0);
    assert_eq!(report.certified_roots(), Vec::<BigInt>::new());
    assert_eq!(report.seeds()[0].seed(), prime.zero());
    assert!(matches!(
        report.seeds()[0].outcome(),
        TorsionXSeedLiftOutcome::SingularModuloPrime
    ));
}

#[test]
fn torsion_x_hensel_lift_reports_simple_seeds_that_do_not_certify_integer_roots() {
    let prime = ReductionPrime::new(5).expect("5 is prime");
    let polynomial =
        IntegerPolynomial::new(vec![BigInt::from(-6), BigInt::from(0), BigInt::from(1)]);
    let report = lift_polynomial_roots_by_hensel(
        3,
        TorsionXPolynomialSource::OddDivisionPolynomial,
        &polynomial,
        prime,
        vec![prime.reduce_u64(1), prime.reduce_u64(4)],
    )
    .expect("uncertified seeds should be non-fatal");

    assert_eq!(report.singular_seed_count(), 0);
    assert_eq!(report.uncertified_seed_count(), 2);
    assert_eq!(report.certified_roots(), Vec::<BigInt>::new());
    assert!(
        report
            .seeds()
            .iter()
            .all(|seed| matches!(seed.outcome(), TorsionXSeedLiftOutcome::NotCertifiedInBound))
    );
}

#[test]
fn reduction_hensel_route_recovers_full_two_torsion_points() {
    let model = integral_model(-1, 0);
    let report = ReducedTorsionLiftReport::from_integral_model(&model)
        .expect("reduction/Hensel route should recover the fixture torsion");

    assert_eq!(report.good_prime().prime().modulus(), 11);
    assert_eq!(
        report.group().shape(),
        RationalTorsionGroupShape::ProductZ2Z2m { m: 1 }
    );
    assert_eq!(
        report.points(),
        &[
            AffinePoint::infinity(),
            rational_point(-1, 0),
            rational_point(0, 0),
            rational_point(1, 0),
        ]
    );
}

#[test]
fn reduction_hensel_route_recovers_cyclic_six_points() {
    let model = integral_model(0, 1);
    let report = ReducedTorsionLiftReport::from_integral_model(&model)
        .expect("reduction/Hensel route should recover cyclic six torsion");

    assert_eq!(
        report.group().shape(),
        RationalTorsionGroupShape::Cyclic { order: 6 }
    );
    assert_eq!(report.points().len(), 6);
    assert!(report.points().contains(&rational_point(2, 3)));
    assert!(report.points().contains(&rational_point(2, -3)));
    for point in report.points() {
        assert!(model.curve().contains(point));
        assert!(
            model
                .curve()
                .exact_mazur_order(point)
                .expect("reported point should be on the curve")
                .is_some()
        );
    }
}

#[test]
fn reduction_hensel_route_can_certify_trivial_torsion_without_lutz_nagell() {
    let model = integral_model(1, 1);
    let report = ReducedTorsionLiftReport::from_integral_model(&model)
        .expect("reduction/Hensel route should finish on trivial torsion fixture");

    assert_eq!(report.group().shape(), RationalTorsionGroupShape::Trivial);
    assert_eq!(report.points(), &[AffinePoint::infinity()]);
}
