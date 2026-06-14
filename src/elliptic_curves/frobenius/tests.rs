use proptest::prelude::*;

use num_bigint::{BigInt, BigUint};
use std::collections::HashSet;
use std::num::NonZeroU32;

use crate::elliptic_curves::frobenius::{
    AbsoluteFrobenius, CharacterSumPointCount, FrobeniusCharacteristicPolynomial,
    FrobeniusCurveType, FrobeniusDiscriminant, FrobeniusLocalZetaFunction,
    FrobeniusTorsionMatrixError, FrobeniusTrace, GroupOrderReport, GroupOrderStrategy,
    HasseGroupOrderStrategy, HasseInterval, MestreConfig, MestreGroupOrderReport, MestreSide,
    MestreStepReport, ModNMatrix2, NTorsionBasis, RelativeFrobenius,
    absolute_frobenius_on_exact_torsion, absolute_frobenius_orbit,
    absolute_frobenius_orbits_on_points, absolute_frobenius_power_point,
    compare_extension_count_with_enumeration, frobenius_matrix_on_n_torsion_basis,
    frobenius_twist_power, relative_frobenius_on_exact_torsion, relative_frobenius_orbit,
    relative_frobenius_orbits_on_points, relative_frobenius_point,
    verify_frobenius_characteristic_equation_at_point,
    verify_frobenius_characteristic_equation_exhaustive, verify_hasse_bound,
    verify_isogeny_frobenius_relation, verify_isogeny_graph_frobenius_relation,
};
use crate::elliptic_curves::{
    AffineCurveModel, AffinePoint, CurveError, CurveModel, EnumerableCurveModel,
    FiniteGroupCurveModel, FrobeniusTraceCurveModel, ShortWeierstrassCurve,
    ShortWeierstrassQuadraticTwist, TwistKind, points_of_exact_order,
};
use crate::fields::{EnumerableFiniteField, Field, FiniteFieldDescriptor, Fp, SqrtField};
use crate::isogenies::graphs::IsogenyGraphBuilder;
use crate::isogenies::{Isogeny, ScalarMultiplicationIsogeny, VeluIsogeny};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::{arb_curve_and_point, arb_nonsingular_curve};
use crate::proptest_support::fields::ProptestF17Sqrt3Field;

type F17 = Fp<17>;
type F43 = Fp<43>;
type F5 = Fp<5>;
type F19 = Fp<19>;
type F41 = Fp<41>;
type F7 = Fp<7>;
type F17Squared = ProptestF17Sqrt3Field;

crate::fields::define_fp_quadratic_extension!(
    spec: F43Sqrt2Spec,
    field: F43Sqrt2,
    base: F43,
    non_residue: 2,
    name: "F43(sqrt(2))",
);

fn alpha() -> <F43Sqrt2 as Field>::Elem {
    F43Sqrt2::element(vec![F43::zero(), F43::one()])
}

fn f41_curve() -> ShortWeierstrassCurve<F41> {
    ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid F7 curve")
}

fn f5_noncyclic_curve() -> ShortWeierstrassCurve<Fp<5>> {
    ShortWeierstrassCurve::<Fp<5>>::new(Fp::<5>::from_i64(-1), Fp::<5>::zero())
        .expect("valid F5 curve")
}

fn depth_one_f41_graph() -> crate::isogenies::graphs::IsogenyGraph<ShortWeierstrassCurve<F41>> {
    IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one F41 graph should build")
}

fn f43_quadratic_twist_factor() -> <F43 as Field>::Elem {
    first_nonsquare::<F43>()
}

fn lift_f17_curve_to_f17_squared(
    curve: &ShortWeierstrassCurve<F17>,
) -> ShortWeierstrassCurve<F17Squared> {
    ShortWeierstrassCurve::<F17Squared>::new(
        F17Squared::from_base(*curve.a()),
        F17Squared::from_base(*curve.b()),
    )
    .expect("lifting a smooth F17 curve to F17^2 should preserve smoothness")
}

fn lift_f43_curve_to_f43_sqrt2(
    curve: &ShortWeierstrassCurve<F43>,
) -> ShortWeierstrassCurve<F43Sqrt2> {
    ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(*curve.a()),
        F43Sqrt2::from_base(*curve.b()),
    )
    .expect("lifting an F43 curve to F43^2 should preserve smoothness")
}

fn find_f43_curve_with_nontrivial_two_torsion_frobenius_basis() -> (
    ShortWeierstrassCurve<F43>,
    ShortWeierstrassCurve<F43Sqrt2>,
    NTorsionBasis<AffinePoint<F43Sqrt2>>,
) {
    let base_curve = ShortWeierstrassCurve::<F43>::new(F43::from_i64(-2), F43::zero())
        .expect("y^2 = x^3 - 2x should be smooth over F43");
    let lifted_curve = lift_f43_curve_to_f43_sqrt2(&base_curve);
    let zero_point = lifted_curve
        .point(F43Sqrt2::zero(), F43Sqrt2::zero())
        .expect("(0,0) should be 2-torsion on the lifted curve");
    let alpha_point = lifted_curve
        .point(alpha(), F43Sqrt2::zero())
        .expect("(sqrt(2),0) should be 2-torsion on the lifted curve");
    let basis = NTorsionBasis::new(&lifted_curve, 2, zero_point, alpha_point)
        .expect("two distinct nonzero 2-torsion points should form a basis");

    (base_curve, lifted_curve, basis)
}

fn non_base_defined_f43_sqrt2_curve_with_two_torsion_basis() -> (
    ShortWeierstrassCurve<F43Sqrt2>,
    NTorsionBasis<AffinePoint<F43Sqrt2>>,
) {
    let u = F43Sqrt2::element(vec![F43::one(), F43::one()]);
    let minus_u_squared = F43Sqrt2::neg(&F43Sqrt2::mul(&u, &u));
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(minus_u_squared, F43Sqrt2::zero())
        .expect("y^2 = x^3 - u^2 x should stay smooth for non-zero u");
    let p0 = curve
        .point(F43Sqrt2::zero(), F43Sqrt2::zero())
        .expect("(0,0) should be 2-torsion");
    let pu = curve
        .point(u.clone(), F43Sqrt2::zero())
        .expect("(u,0) should be 2-torsion");

    let basis = NTorsionBasis::new(&curve, 2, p0, pu)
        .expect("two distinct nonzero 2-torsion points should form a basis");

    (curve, basis)
}

fn nz(n: u32) -> core::num::NonZeroU32 {
    core::num::NonZeroU32::new(n).expect("test degrees are positive")
}

fn find_non_fixed_point(curve: &ShortWeierstrassCurve<F43Sqrt2>) -> AffinePoint<F43Sqrt2> {
    for x in F43Sqrt2::elements() {
        for y in F43Sqrt2::elements() {
            if let Ok(point) = curve.point(x.clone(), y) {
                let image = absolute_frobenius_power_point(curve, &point, 1)
                    .expect("absolute Frobenius should evaluate on on-curve inputs");
                if image != point {
                    return point;
                }
            }
        }
    }

    panic!("expected a point in E(F43^2) that is not fixed by Frob_43");
}

fn first_nonsquare<F>() -> F::Elem
where
    F: EnumerableFiniteField + SqrtField,
{
    F::elements()
        .into_iter()
        .find(|value| !F::is_zero(value) && !F::has_square_root(value))
        .expect("small odd prime fields should contain non-squares")
}

#[test]
fn absolute_frobenius_metadata_uses_the_prime_characteristic() {
    let frobenius = AbsoluteFrobenius::for_field::<F43>(3);

    assert_eq!(frobenius.characteristic, 43);
    assert_eq!(frobenius.power, 3);
    assert!(!frobenius.is_identity());
}

#[test]
fn relative_frobenius_metadata_uses_the_full_base_field_descriptor() {
    let frobenius = RelativeFrobenius::for_field::<F43Sqrt2>(2);

    assert_eq!(frobenius.base_field.characteristic, 43);
    assert_eq!(frobenius.base_field.extension_degree.get(), 2);
    assert_eq!(frobenius.power, 2);
}

#[test]
fn frobenius_twist_can_change_extension_field_coefficients() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(alpha(), F43Sqrt2::zero())
        .expect("non-zero linear coefficient should define a smooth curve");

    let first_twist = frobenius_twist_power(&curve, 1).expect("Frobenius twist should be valid");
    let second_twist =
        frobenius_twist_power(&curve, 2).expect("second Frobenius twist should be valid");

    assert_ne!(first_twist, curve);
    assert_eq!(second_twist, curve);
}

#[test]
fn absolute_frobenius_fixes_the_point_at_infinity() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = AffinePoint::<F43>::infinity();

    let image = absolute_frobenius_power_point(&curve, &point, 4)
        .expect("Frobenius should fix the distinguished identity");

    assert_eq!(image, point);
}

#[test]
fn absolute_frobenius_fixes_prime_field_rational_points() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) should lie on y^2 = x^3 + x + 1 over F43");

    let image = absolute_frobenius_power_point(&curve, &point, 1)
        .expect("Frobenius should evaluate on rational points");

    assert_eq!(image, point);
}

#[test]
fn relative_frobenius_orbit_on_a_rational_point_is_a_singleton() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) should lie on the curve");

    let orbit = relative_frobenius_orbit(&curve, &point).expect("relative orbit should evaluate");

    assert_eq!(orbit.start(), &point);
    assert_eq!(orbit.points(), &[point]);
    assert_eq!(orbit.period(), 1);
}

#[test]
fn relative_frobenius_orbits_on_points_are_singletons() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let orbits =
        relative_frobenius_orbits_on_points(&curve).expect("relative point orbits should evaluate");

    assert_eq!(orbits.len(), curve.order());
    assert!(orbits.iter().all(|orbit| orbit.period() == 1));
    assert_eq!(
        orbits.iter().map(|orbit| orbit.period()).sum::<usize>(),
        curve.order()
    );
}

#[test]
fn relative_frobenius_can_fix_a_point_that_absolute_frobenius_moves() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let point = find_non_fixed_point(&curve);

    let absolute_image = absolute_frobenius_power_point(&curve, &point, 1)
        .expect("absolute Frobenius should evaluate");
    let relative_image =
        relative_frobenius_point(&curve, &point).expect("relative Frobenius should evaluate");

    assert_ne!(absolute_image, point);
    assert_eq!(relative_image, point);
}

#[test]
fn rational_two_torsion_basis_over_the_base_field_gives_the_identity_matrix() {
    let curve = f5_noncyclic_curve();
    let two_torsion =
        points_of_exact_order(&curve, 2).expect("F5 curve should have full rational 2-torsion");
    let basis = NTorsionBasis::new(&curve, 2, two_torsion[0].clone(), two_torsion[1].clone())
        .expect("two independent rational 2-torsion points should form a basis");
    let trace = curve
        .frobenius_trace()
        .expect("small enumerable F5 curve should supply a Frobenius trace");

    let report = frobenius_matrix_on_n_torsion_basis(&curve, trace, basis)
        .expect("matrix report should build on rational 2-torsion");

    assert_eq!(report.matrix().entries(), [[1, 0], [0, 1]]);
    assert!(report.trace_matches_mod_n());
    assert!(report.determinant_matches_mod_n());
}

#[test]
fn nontrivial_extension_two_torsion_basis_still_matches_trace_and_degree_mod_n() {
    let (base_curve, lifted_curve, basis) =
        find_f43_curve_with_nontrivial_two_torsion_frobenius_basis();
    let trace = base_curve
        .frobenius_trace()
        .expect("base F43 curve should supply a Frobenius trace");

    let report = frobenius_matrix_on_n_torsion_basis(&lifted_curve, trace.clone(), basis)
        .expect("matrix report should build over the lifted curve");

    assert_ne!(report.matrix().entries(), [[1, 0], [0, 1]]);
    assert_eq!(report.matrix().modulus(), 2);
    assert!(report.trace_matches_mod_n());
    assert!(report.determinant_matches_mod_n());
    assert_eq!(report.frobenius_trace(), &trace);
}

#[test]
fn torsion_basis_rejects_dependent_points() {
    let curve = f5_noncyclic_curve();
    let two_torsion =
        points_of_exact_order(&curve, 2).expect("F5 curve should have full rational 2-torsion");

    assert_eq!(
        NTorsionBasis::new(&curve, 2, two_torsion[0].clone(), two_torsion[0].clone()),
        Err(FrobeniusTorsionMatrixError::DependentTorsionBasis)
    );
}

#[test]
fn torsion_basis_rejects_orders_divisible_by_the_characteristic() {
    let curve = f5_noncyclic_curve();
    let point = curve.identity();

    assert_eq!(
        NTorsionBasis::new(&curve, 5, point.clone(), point),
        Err(
            FrobeniusTorsionMatrixError::CharacteristicDividesTorsionOrder {
                characteristic: 5,
                order: 5,
            }
        )
    );
}

#[test]
fn matrix_report_rejects_traces_over_the_wrong_characteristic() {
    let curve = f5_noncyclic_curve();
    let two_torsion =
        points_of_exact_order(&curve, 2).expect("F5 curve should have full rational 2-torsion");
    let basis = NTorsionBasis::new(&curve, 2, two_torsion[0].clone(), two_torsion[1].clone())
        .expect("two independent rational 2-torsion points should form a basis");
    let wrong_trace =
        FrobeniusTrace::from_order(FiniteFieldDescriptor::new(43, nz(1)).unwrap(), 41)
            .expect("t = 3 should be valid over F43");

    assert_eq!(
        frobenius_matrix_on_n_torsion_basis(&curve, wrong_trace, basis),
        Err(
            FrobeniusTorsionMatrixError::TraceBaseFieldCharacteristicMismatch {
                trace_characteristic: 43,
                curve_characteristic: 5,
            }
        )
    );
}

#[test]
fn matrix_report_rejects_curves_not_fixed_by_the_trace_frobenius_power() {
    let (curve, basis) = non_base_defined_f43_sqrt2_curve_with_two_torsion_basis();
    let trace = FrobeniusTrace::from_order(FiniteFieldDescriptor::new(43, nz(1)).unwrap(), 41)
        .expect("t = 3 should be valid over F43");

    assert_eq!(
        frobenius_matrix_on_n_torsion_basis(&curve, trace, basis),
        Err(
            FrobeniusTorsionMatrixError::FrobeniusTraceDoesNotPreserveCurve {
                extension_degree: 1,
            }
        )
    );
}

#[test]
fn matrix_helper_reduces_trace_and_determinant_mod_n() {
    let matrix = ModNMatrix2::new(5, [[4, 3], [2, 1]]).expect("entries should already be reduced");

    assert_eq!(matrix.trace_mod_n(), 0);
    assert_eq!(matrix.determinant_mod_n(), Ok(3));
}

#[test]
fn absolute_frobenius_orbit_can_have_period_two_over_f43_squared() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let point = find_non_fixed_point(&curve);

    let orbit =
        absolute_frobenius_orbit(&curve, &point, 1).expect("absolute orbit should evaluate");
    let image =
        absolute_frobenius_power_point(&curve, &point, 1).expect("absolute Frobenius should apply");

    assert_eq!(orbit.start(), &point);
    assert_eq!(orbit.points(), &[point.clone(), image]);
    assert_eq!(orbit.period(), 2);
}

#[test]
fn absolute_frobenius_orbits_on_points_partition_a_quadratic_extension_curve() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let orbits = absolute_frobenius_orbits_on_points(&curve, 1)
        .expect("absolute Frobenius point orbits should evaluate");

    assert_eq!(
        orbits.iter().map(|orbit| orbit.period()).sum::<usize>(),
        curve.order()
    );
    assert!(orbits.iter().any(|orbit| orbit.period() == 2));
}

#[test]
fn absolute_frobenius_orbit_rejects_curves_not_preserved_by_the_chosen_power() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(alpha(), F43Sqrt2::zero())
        .expect("valid extension-defined curve");
    let point = AffinePoint::<F43Sqrt2>::infinity();

    assert_eq!(
        absolute_frobenius_orbit(&curve, &point, 1),
        Err(crate::elliptic_curves::CurveError::AbsoluteFrobeniusDoesNotPreserveCurve { power: 1 })
    );
}

#[test]
fn curves_over_enumerable_quadratic_extensions_can_list_points_and_orders() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let points = curve.points();

    assert!(!points.is_empty());
    assert_eq!(curve.order(), points.len());
}

#[test]
fn frobenius_trace_matches_the_counting_formula() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let report = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    assert_eq!(report.base_field().characteristic, 43);
    assert_eq!(report.base_field().extension_degree.get(), 1);
    assert_eq!(report.field_order(), 43);
    assert_eq!(report.curve_order(), curve.order() as u64);
    assert_eq!(report.trace(), 43_i64 + 1 - curve.order() as i64);
}

#[test]
fn frobenius_trace_from_order_and_order_from_trace_roundtrip() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let report = FrobeniusTrace::from_order(base_field.clone(), 48)
        .expect("small Frobenius trace package should build");

    assert_eq!(report.trace(), -4);
    assert_eq!(
        FrobeniusTrace::curve_order_from_trace(base_field, report.trace()),
        Ok(report.curve_order())
    );
}

#[test]
fn frobenius_curve_order_from_trace_rejects_non_positive_orders() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    assert_eq!(
        FrobeniusTrace::curve_order_from_trace(base_field, 100),
        Err(crate::elliptic_curves::CurveError::InvalidFrobeniusTrace { trace: 100 })
    );
}

#[test]
fn extension_count_over_degree_one_recovers_the_base_curve_order() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let report = trace.curve_order_over_extension(nz(1));

    assert_eq!(report.frobenius_trace(), &trace);
    assert_eq!(report.extension_degree(), nz(1));
    assert_eq!(report.extension_field_order(), &BigUint::from(43u32));
    assert_eq!(report.power_sum(), &BigInt::from(trace.trace()));
    assert_eq!(report.curve_order(), &BigUint::from(curve.order() as u64));
}

#[test]
fn extension_count_over_degree_two_matches_explicit_quadratic_extension_enumeration() {
    let base_curve =
        ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let extension_curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let trace = base_curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let report = trace.curve_order_over_extension(nz(2));
    let expected_power_sum =
        BigInt::from(trace.trace()) * BigInt::from(trace.trace()) - BigInt::from(2 * 43);

    assert_eq!(report.extension_degree(), nz(2));
    assert_eq!(report.extension_field_order(), &BigUint::from(43u32).pow(2));
    assert_eq!(report.power_sum(), &expected_power_sum);
    assert_eq!(
        report.curve_order(),
        &BigUint::from(extension_curve.order() as u64)
    );
}

#[test]
fn extension_count_sequence_matches_individual_reports() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let sequence = trace.curve_orders_over_extensions_through(nz(4));

    assert_eq!(sequence.frobenius_trace(), &trace);
    assert_eq!(sequence.reports().len(), 4);
    assert_eq!(
        sequence
            .reports()
            .iter()
            .map(|report| report.extension_degree().get())
            .collect::<Vec<_>>(),
        vec![1, 2, 3, 4]
    );

    for report in sequence.reports() {
        let individual = trace.curve_order_over_extension(report.extension_degree());
        assert_eq!(report, &individual);
    }
}

#[test]
fn extension_count_can_represent_degrees_far_beyond_u128_field_powers() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let report = trace.curve_order_over_extension(nz(25));

    assert!(report.extension_field_order() > &BigUint::from(u128::MAX));
    assert!(report.curve_order() > &BigUint::from(u128::MAX));
}

#[test]
fn extension_count_comparison_report_distinguishes_frobenius_and_enumeration_routes() {
    let base_curve =
        ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid base curve");
    let extension_curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::one()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");
    let trace = base_curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let comparison = compare_extension_count_with_enumeration(&extension_curve, &trace)
        .expect("comparison between Frobenius and exhaustive paths should compute");

    assert_eq!(comparison.trace_base_field(), trace.base_field());
    assert_eq!(comparison.curve_base_field().characteristic, 43);
    assert_eq!(comparison.curve_base_field().extension_degree.get(), 2);
    assert_eq!(comparison.relative_extension_degree(), nz(2));
    assert_eq!(
        comparison.frobenius_count().curve_order(),
        comparison.exhaustive_curve_order()
    );
    assert_eq!(
        comparison.exhaustive_curve_order(),
        &BigUint::from(extension_curve.order() as u64)
    );
    assert!(comparison.agrees());
}

#[test]
fn characteristic_polynomial_is_derived_from_the_frobenius_trace() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let polynomial = trace.characteristic_polynomial();

    assert_eq!(polynomial.base_field(), trace.base_field());
    assert_eq!(polynomial.trace(), trace.trace());
    assert_eq!(polynomial.field_order(), trace.field_order());
}

#[test]
fn characteristic_polynomial_discriminant_and_evaluation_match_the_formula() {
    let polynomial = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        -4,
    );

    assert_eq!(polynomial.discriminant(), 16 - 4 * 43);
    assert_eq!(polynomial.evaluate_at_integer(0), 43);
    assert_eq!(polynomial.evaluate_at_integer(1), 48);
    assert_eq!(polynomial.evaluate_at_integer(2), 55);
}

#[test]
fn frobenius_discriminant_is_derived_from_the_trace_package() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let discriminant = trace.discriminant();

    assert_eq!(discriminant.frobenius_trace(), &trace);
    assert_eq!(discriminant.base_field(), trace.base_field());
    assert_eq!(discriminant.curve_order(), trace.curve_order());
    assert_eq!(discriminant.trace(), trace.trace());
}

#[test]
fn frobenius_discriminant_matches_the_quadratic_formula() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace = FrobeniusTrace::from_order(base_field, 41).expect("t = 3 should be valid over F43");

    let discriminant = trace.discriminant();

    assert_eq!(
        discriminant.quadratic_discriminant().value(),
        &num_bigint::BigInt::from(-163)
    );
    assert!(discriminant.is_negative());
    assert!(!discriminant.is_zero());
    assert!(!discriminant.is_positive());
    assert!(discriminant.is_fundamental());
}

#[test]
fn frobenius_discriminant_constructor_matches_trace_method() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");

    let from_method = trace.discriminant();
    let from_constructor = FrobeniusDiscriminant::new(trace.clone());

    assert_eq!(from_constructor, from_method);
    assert_eq!(
        from_constructor.quadratic_discriminant().value(),
        &num_bigint::BigInt::from(-156)
    );
}

#[test]
fn frobenius_discriminant_quadratic_factorization_matches_the_expected_split() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");
    let discriminant = trace.discriminant();

    let factorization = discriminant.quadratic_factorization().expect(
        "negative Frobenius discriminant should admit an imaginary quadratic factorization",
    );

    assert_eq!(
        factorization.discriminant(),
        discriminant.quadratic_discriminant()
    );
    assert_eq!(factorization.conductor(), &num_bigint::BigUint::from(2u8));
    assert_eq!(
        factorization.fundamental_discriminant(),
        &crate::elliptic_curves::QuadraticDiscriminant::new(-39)
    );
}

#[test]
fn frobenius_discriminant_can_build_the_imaginary_quadratic_order() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");
    let discriminant = trace.discriminant();

    let order = discriminant
        .frobenius_order()
        .expect("negative Frobenius discriminant should define an imaginary quadratic order");

    assert_eq!(
        order.fundamental_discriminant(),
        &crate::elliptic_curves::QuadraticDiscriminant::new(-39)
    );
    assert_eq!(order.conductor(), &num_bigint::BigUint::from(2u8));
    assert_eq!(order.discriminant(), discriminant.quadratic_discriminant());
    assert!(!order.is_maximal());
}

#[test]
fn frobenius_order_matches_the_existing_imaginary_quadratic_order_helper() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");
    let discriminant = trace.discriminant();

    let from_existing = discriminant
        .frobenius_order()
        .expect("existing helper should build the order");
    let from_frobenius_order = discriminant
        .frobenius_order()
        .expect("frobenius_order should build the same order");

    assert_eq!(from_frobenius_order, from_existing);
}

#[test]
fn frobenius_order_is_contained_in_the_maximal_order() {
    let base_field =
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid");
    let trace =
        FrobeniusTrace::from_order(base_field, 48).expect("t = -4 should be valid over F43");
    let discriminant = trace.discriminant();
    let frobenius_order = discriminant
        .frobenius_order()
        .expect("Frobenius order should exist in the imaginary case");
    let maximal_order = discriminant
        .maximal_order()
        .expect("the same factorization should recover O_K");

    assert!(frobenius_order.is_suborder_of(&maximal_order));
    assert!(maximal_order.is_overorder_of(&frobenius_order));
}

#[test]
fn endomorphism_ring_report_from_frobenius_is_honest_about_candidates() {
    let base_field = FiniteFieldDescriptor::new(13, nz(1))
        .expect("F_13 metadata should be internally consistent");
    let trace =
        FrobeniusTrace::from_order(base_field, 18).expect("t = -4 should be valid over F13");
    let report = trace
        .discriminant()
        .endomorphism_ring_report()
        .expect("ordinary Frobenius data should produce a candidate report");

    assert_eq!(report.frobenius_discriminant().trace(), -4);
    assert!(matches!(
        report,
        crate::elliptic_curves::EndomorphismRingReport::OrdinaryQuadraticOrderCandidates { .. }
    ));
    assert_eq!(report.candidate_count(), Some(2));
    assert_eq!(report.sandwich_inclusion_holds(), Some(true));
    assert_eq!(
        report
            .candidate_orders()
            .expect("ordinary branch should expose candidate orders")
            .iter()
            .map(|order| order.conductor().clone())
            .collect::<Vec<_>>(),
        vec![BigUint::from(1u8), BigUint::from(3u8)]
    );
}

#[test]
fn supersingular_frobenius_data_uses_the_supersingular_report_branch() {
    let base_field =
        FiniteFieldDescriptor::new(5, nz(1)).expect("F_5 metadata should be internally consistent");
    let trace = FrobeniusTrace::from_order(base_field, 6).expect("t = 0 should be valid over F5");
    let report = trace
        .discriminant()
        .endomorphism_ring_report()
        .expect("supersingular Frobenius data should still produce a report");

    assert!(matches!(
        report,
        crate::elliptic_curves::EndomorphismRingReport::SupersingularQuaternionicPlaceholder { .. }
    ));
    assert_eq!(report.frobenius_discriminant().trace(), 0);
    assert_eq!(report.candidate_count(), None);
}

#[test]
fn characteristic_polynomial_pretty_and_display_use_the_same_educational_surface() {
    let positive_trace = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        3,
    );
    let zero_trace = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        0,
    );
    let negative_trace = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        -4,
    );

    assert_eq!(positive_trace.pretty(), "T^2 - 3T + 43");
    assert_eq!(zero_trace.pretty(), "T^2 + 43");
    assert_eq!(negative_trace.pretty(), "T^2 + 4T + 43");
    assert_eq!(format!("{negative_trace}"), negative_trace.pretty());
}

#[test]
fn local_zeta_function_is_derived_from_the_characteristic_polynomial() {
    let polynomial = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F43 descriptor should be valid"),
        -4,
    );

    let zeta = polynomial.local_zeta_function();

    assert_eq!(zeta.characteristic_polynomial(), &polynomial);
    assert_eq!(zeta.base_field(), polynomial.base_field());
    assert_eq!(zeta.field_order(), polynomial.field_order());
    assert_eq!(zeta.trace(), polynomial.trace());
}

#[test]
fn local_zeta_function_is_derived_from_the_frobenius_trace() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let zeta = trace.local_zeta_function();

    assert_eq!(
        zeta,
        trace.characteristic_polynomial().local_zeta_function()
    );
}

#[test]
fn local_zeta_function_pretty_and_display_use_the_same_educational_surface() {
    let prime_field_zeta = FrobeniusLocalZetaFunction::from_characteristic_polynomial(
        FrobeniusCharacteristicPolynomial::new(
            FiniteFieldDescriptor::new(17, core::num::NonZeroU32::new(1).expect("1 is non-zero"))
                .expect("F17 descriptor should be valid"),
            2,
        ),
    );
    let extension_field_zeta = FrobeniusLocalZetaFunction::from_characteristic_polynomial(
        FrobeniusCharacteristicPolynomial::new(
            FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(2).expect("2 is non-zero"))
                .expect("F43^2 descriptor should be valid"),
            -4,
        ),
    );

    assert_eq!(prime_field_zeta.numerator_pretty(), "1 - 2T + 17T²");
    assert_eq!(prime_field_zeta.denominator_pretty(), "(1 - T)(1 - 17T)");
    assert_eq!(
        prime_field_zeta.pretty(),
        "Z(E/F_17, T) = (1 - 2T + 17T²) / ((1 - T)(1 - 17T))"
    );
    assert_eq!(format!("{prime_field_zeta}"), prime_field_zeta.pretty());

    assert_eq!(
        extension_field_zeta.pretty(),
        "Z(E/F_(43^2), T) = (1 + 4T + 1849T²) / ((1 - T)(1 - 1849T))"
    );
}

#[test]
fn scalar_isogeny_frobenius_relation_holds_on_a_small_f41_curve() {
    let isogeny = ScalarMultiplicationIsogeny::new(f41_curve(), 2)
        .expect("scalar multiplication isogeny should build");

    let relation = verify_isogeny_frobenius_relation(&isogeny)
        .expect("Frobenius relation should compute on both sides of the isogeny");

    assert_eq!(relation.degree(), isogeny.degree());
    assert!(relation.same_curve_order());
    assert!(relation.same_trace());
    assert!(relation.holds());
}

#[test]
fn velu_isogeny_frobenius_relation_holds_on_the_small_f41_example() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample generator should lie on the curve");
    let isogeny =
        VeluIsogeny::from_generator(domain, generator).expect("sample Vélu isogeny should build");

    let relation = verify_isogeny_frobenius_relation(&isogeny)
        .expect("Frobenius relation should compute on the Vélu example");

    assert_eq!(relation.degree(), isogeny.degree());
    assert_eq!(
        relation.domain().curve_order(),
        relation.codomain().curve_order()
    );
    assert_eq!(relation.domain().trace(), relation.codomain().trace());
    assert!(relation.same_curve_order());
    assert!(relation.same_trace());
    assert!(relation.holds());
}

#[test]
fn isogeny_graph_frobenius_report_holds_on_the_depth_one_f41_graph() {
    let graph = depth_one_f41_graph();

    let report = verify_isogeny_graph_frobenius_relation(&graph)
        .expect("graph Frobenius relation should compute on the F41 graph");

    assert_eq!(report.checked_nodes(), graph.node_count());
    assert_eq!(report.checked_edges(), graph.edge_count());
    assert_eq!(report.reference_node().0, 0);
    assert!(report.all_same_curve_order());
    assert!(report.all_same_trace());
    assert!(report.holds());
    assert!(report.node_reports().iter().all(|node| node.holds()));
}

#[test]
fn isogeny_graph_frobenius_report_handles_a_singleton_graph() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("singleton F41 graph should build");

    let report = verify_isogeny_graph_frobenius_relation(&graph)
        .expect("singleton graph Frobenius relation should compute");

    assert_eq!(report.checked_nodes(), 1);
    assert_eq!(report.checked_edges(), 0);
    assert!(report.all_same_curve_order());
    assert!(report.all_same_trace());
    assert!(report.holds());
    assert_eq!(report.node_reports().len(), 1);
    assert!(report.node_reports()[0].holds());
}

#[test]
fn quadratic_twist_frobenius_relation_holds_for_a_nontrivial_f19_twist() {
    let curve = ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3))
        .expect("valid F19 curve");
    let package = ShortWeierstrassQuadraticTwist::new(curve, first_nonsquare::<F19>())
        .expect("non-zero twist factor should define a quadratic-twist package");

    let relation = package
        .frobenius_relation()
        .expect("Frobenius relation should compute by exhaustive counting");

    assert_eq!(relation.sum_orders(), 2 * 19 + 2);
    assert_eq!(relation.expected_sum(), 2 * 19 + 2);
    assert_eq!(relation.twist_kind(), TwistKind::Quadratic);
    assert!(relation.holds());
    assert!(relation.trace_negation_holds());
    assert!(relation.matches_twist_kind_expectation());
}

#[test]
fn quadratic_twist_frobenius_relation_reports_trivial_twists_without_erroring() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid F43 curve");
    let package = ShortWeierstrassQuadraticTwist::new(curve, F43::from_i64(4))
        .expect("square twist factor should still define a twist package");

    let relation = package
        .frobenius_relation()
        .expect("trivial twists should still produce a report");

    assert_eq!(relation.twist_kind(), TwistKind::Trivial);
    assert_eq!(
        relation.original().curve_order(),
        relation.twist().curve_order()
    );
    assert_eq!(relation.original().trace(), relation.twist().trace());
    assert!(relation.same_curve_order_holds());
    assert!(relation.same_trace_holds());
    assert!(!relation.trace_negation_holds());
    assert!(!relation.holds());
    assert!(relation.matches_twist_kind_expectation());
}

#[test]
fn j_1728_f43_nonsquare_factor_can_still_be_trivial_and_frobenius_equal() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::zero())
        .expect("y^2 = x^3 + x should be smooth over F43");
    let factor = F43::from_i64(2);
    let package = ShortWeierstrassQuadraticTwist::new(curve, factor)
        .expect("non-zero factor should define a twist package");
    let relation = package
        .frobenius_relation()
        .expect("Frobenius relation should compute on the exceptional j = 1728 case");

    assert!(!F43::has_square_root(package.factor()));
    assert_eq!(package.kind(), TwistKind::Trivial);
    assert_eq!(relation.twist_kind(), TwistKind::Trivial);
    assert_eq!(relation.original().curve_order(), 44);
    assert_eq!(relation.twist().curve_order(), 44);
    assert_eq!(relation.original().trace(), 0);
    assert_eq!(relation.twist().trace(), 0);
    assert!(relation.same_curve_order_holds());
    assert!(relation.same_trace_holds());
    assert!(relation.trace_negation_holds());
    assert!(relation.holds());
    assert!(relation.matches_twist_kind_expectation());
}

#[test]
fn hasse_bound_report_matches_the_frobenius_trace_data() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let expected_trace = 43_i64 + 1 - curve.order() as i64;
    let expected_trace_square = i128::from(expected_trace) * i128::from(expected_trace);

    let report = verify_hasse_bound(&curve).expect("Hasse bound should verify");

    assert_eq!(report.frobenius_trace().field_order(), 43);
    assert_eq!(report.frobenius_trace().curve_order(), curve.order() as u64);
    assert_eq!(report.frobenius_trace().trace(), expected_trace);
    assert_eq!(report.trace_square(), expected_trace_square);
    assert_eq!(report.bound_square(), 4 * 43);
    assert_eq!(report.slack(), 4 * 43 - expected_trace_square);
    assert!(report.holds());
}

#[test]
fn hasse_interval_for_non_square_q_matches_the_discrete_formula() {
    let interval = HasseInterval::for_q(43).expect("q = 43 should define H(q)");

    assert_eq!(interval.q(), 43);
    assert_eq!(interval.lower(), 31);
    assert_eq!(interval.upper(), 57);
    assert_eq!(interval.span(), 26);
    assert_eq!(interval.candidate_count(), 27);
    assert!(interval.contains(31));
    assert!(interval.contains(41));
    assert!(interval.contains(57));
    assert!(!interval.contains(30));
    assert!(!interval.contains(58));
    assert_eq!(interval.as_range_inclusive().collect::<Vec<_>>()[0], 31);
}

#[test]
fn hasse_interval_multiple_helpers_distinguish_none_unique_and_several() {
    let interval = HasseInterval::for_q(43).expect("q = 43 should define H(q)");

    assert_eq!(interval.first_multiple_of(11), Some(33));
    assert_eq!(interval.last_multiple_of(11), Some(55));
    assert_eq!(interval.multiple_count_of(11), 3);
    assert_eq!(interval.multiples_of(11), vec![33, 44, 55]);
    assert_eq!(interval.unique_multiple_of(11), None);

    assert_eq!(interval.first_multiple_of(23), Some(46));
    assert_eq!(interval.last_multiple_of(23), Some(46));
    assert_eq!(interval.multiple_count_of(23), 1);
    assert_eq!(interval.unique_multiple_of(23), Some(46));

    assert_eq!(interval.first_multiple_of(29), None);
    assert_eq!(interval.last_multiple_of(29), None);
    assert_eq!(interval.multiple_count_of(29), 0);
    assert_eq!(interval.multiples_of(29), Vec::<u128>::new());
    assert_eq!(interval.unique_multiple_of(29), None);

    assert_eq!(interval.first_multiple_of(0), None);
    assert_eq!(interval.last_multiple_of(0), None);
    assert_eq!(interval.multiple_count_of(0), 0);
    assert_eq!(interval.multiples_of(0), Vec::<u128>::new());
    assert_eq!(interval.unique_multiple_of(0), None);
}

#[test]
fn naive_hasse_multiple_search_finds_the_first_annihilating_multiple_in_h7() {
    let curve = f7_curve();
    let point = curve
        .generator()
        .expect("the F7 sample curve should be cyclic");

    let report = curve
        .find_annihilating_multiple_in_hasse_interval_naive(&point)
        .expect("naive Hasse search should succeed on on-curve inputs");

    assert_eq!(report.q(), 7);
    assert_eq!(
        report.interval(),
        &HasseInterval::for_q(7).expect("H(7) should exist")
    );
    assert_eq!(report.first_annihilating_multiple(), Some(6));
    assert_eq!(report.tested_candidates(), 4);
    assert_eq!(
        report.steps().first().map(|step| step.candidate_multiple()),
        Some(3)
    );
    assert_eq!(
        report.steps().last().map(|step| step.candidate_multiple()),
        Some(6)
    );
    assert!(
        curve.is_identity(
            report
                .steps()
                .last()
                .expect("the first annihilating candidate should be recorded")
                .image()
        )
    );
}

#[test]
fn naive_hasse_multiple_search_kills_the_identity_at_the_lower_endpoint() {
    let curve = f41_curve();
    let identity = curve.identity();

    let report = curve
        .find_annihilating_multiple_in_hasse_interval_naive(&identity)
        .expect("the identity should be searchable");

    assert_eq!(report.first_annihilating_multiple(), Some(30));
    assert_eq!(report.tested_candidates(), 1);
    assert_eq!(report.steps().len(), 1);
    assert!(curve.is_identity(report.steps()[0].image()));
}

#[test]
fn naive_hasse_multiple_search_rejects_off_curve_inputs() {
    let curve = f41_curve();
    let invalid = AffinePoint::<F41>::new(F41::from_i64(2), F41::from_i64(2));

    assert_eq!(
        curve.find_annihilating_multiple_in_hasse_interval_naive(&invalid),
        Err(crate::elliptic_curves::CurveError::PointNotOnCurve)
    );
}

#[test]
fn hasse_interval_from_trace_matches_trace_method() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");

    let from_trace = HasseInterval::from_trace(&trace);
    let from_method = trace.hasse_interval();

    assert_eq!(from_trace, from_method);
    assert!(from_trace.contains(trace.curve_order().into()));
}

#[test]
fn hasse_interval_rejects_invalid_field_orders() {
    assert_eq!(
        HasseInterval::for_q(0),
        Err(crate::elliptic_curves::CurveError::InvalidHasseIntervalFieldOrder { field_order: 0 })
    );
    assert_eq!(
        HasseInterval::for_q(1),
        Err(crate::elliptic_curves::CurveError::InvalidHasseIntervalFieldOrder { field_order: 1 })
    );
}

#[test]
fn character_sum_count_matches_exhaustive_order_and_trace_over_f43() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let report = curve
        .group_order_by_quadratic_character()
        .expect("character-sum count should succeed");
    let trace = curve
        .frobenius_trace()
        .expect("exhaustive Frobenius trace should compute");

    assert_eq!(report.base_field(), trace.base_field());
    assert_eq!(report.field_order(), 43);
    assert_eq!(report.curve_order(), curve.order() as u128);
    assert_eq!(report.trace(), i128::from(trace.trace()));
    assert_eq!(
        report.curve_order() as i128,
        report.field_order() as i128 + 1 + report.character_sum()
    );
    assert!(report.hasse_interval().contains(report.curve_order()));
    assert_eq!(
        report
            .to_frobenius_trace()
            .expect("report should convert to the shared Frobenius trace"),
        trace
    );
}

#[test]
fn character_sum_count_matches_exhaustive_order_and_trace_over_quadratic_extension() {
    let base_curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
        .expect("valid F17 curve");
    let curve = lift_f17_curve_to_f17_squared(&base_curve);

    let report = curve
        .group_order_by_quadratic_character()
        .expect("character-sum count should succeed");
    let trace = curve
        .frobenius_trace()
        .expect("exhaustive Frobenius trace should compute");

    assert_eq!(report.base_field(), trace.base_field());
    assert_eq!(report.field_order(), 17_u128.pow(2));
    assert_eq!(report.curve_order(), curve.order() as u128);
    assert_eq!(report.trace(), i128::from(trace.trace()));
    assert_eq!(
        report
            .to_frobenius_trace()
            .expect("report should convert to the shared Frobenius trace"),
        trace
    );
}

#[test]
fn character_sum_report_constructor_recovers_trace_by_negating_the_sum() {
    let base_field = FiniteFieldDescriptor::new(43, nz(1)).expect("F43 metadata should be valid");
    let report =
        CharacterSumPointCount::new(base_field.clone(), -3).expect("character-sum report builds");

    assert_eq!(report.base_field(), &base_field);
    assert_eq!(report.field_order(), 43);
    assert_eq!(report.character_sum(), -3);
    assert_eq!(report.curve_order(), 41);
    assert_eq!(report.trace(), 3);
}

#[test]
fn unified_group_order_api_exposes_both_routes() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let automatic = curve
        .group_order_by(GroupOrderStrategy::Auto)
        .expect("automatic group order should succeed");
    let exhaustive = curve
        .group_order_by(GroupOrderStrategy::Exhaustive)
        .expect("exhaustive group order should succeed");

    assert_eq!(automatic.strategy(), GroupOrderStrategy::QuadraticCharacter);
    assert_eq!(exhaustive.strategy(), GroupOrderStrategy::Exhaustive);
    assert_eq!(automatic.curve_order(), exhaustive.curve_order());
    assert_eq!(automatic.trace(), exhaustive.trace());
    assert_eq!(
        automatic
            .to_frobenius_trace()
            .expect("automatic report should convert to a shared trace"),
        curve
            .frobenius_trace_by(GroupOrderStrategy::Exhaustive)
            .expect("exhaustive trace should compute")
    );
}

#[test]
fn unified_group_order_report_preserves_the_underlying_variant() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let report = curve
        .group_order_by(GroupOrderStrategy::QuadraticCharacter)
        .expect("quadratic-character group order should succeed");

    match report {
        GroupOrderReport::QuadraticCharacter(character_sum) => {
            assert_eq!(character_sum.curve_order(), curve.order() as u128);
        }
        GroupOrderReport::ExhaustiveTrace(_) => {
            panic!("quadratic-character count should preserve its variant")
        }
        GroupOrderReport::MestreFp(_) => {
            panic!("quadratic-character count should not use the Mestre variant")
        }
        GroupOrderReport::FromExponentLowerBound(_) => {
            panic!("quadratic-character count should not use the lower-bound variant")
        }
    }
}

#[test]
fn exponent_lower_bound_route_can_recover_one_unique_group_order() {
    let curve = ShortWeierstrassCurve::<F5>::new(F5::zero(), F5::one()).expect("valid curve");

    let report = curve
        .group_order_by(GroupOrderStrategy::FromExponentLowerBoundAndPointCount {
            exponent_lower_bound: BigUint::from(6u8),
            hasse_strategy: HasseGroupOrderStrategy::Exhaustive,
        })
        .expect("unique H(q) multiple should recover the group order");

    assert_eq!(
        report.strategy(),
        GroupOrderStrategy::FromExponentLowerBoundAndPointCount {
            exponent_lower_bound: BigUint::from(6u8),
            hasse_strategy: HasseGroupOrderStrategy::Exhaustive,
        }
    );
    assert_eq!(report.curve_order(), 6);

    let GroupOrderReport::FromExponentLowerBound(verification) = report else {
        panic!("expected the lower-bound route to preserve its own report variant");
    };

    assert_eq!(verification.exponent_lower_bound(), &BigUint::from(6u8));
    assert_eq!(verification.verified_group_order(), Some(6));
}

#[test]
fn mestre_group_order_report_preserves_original_curve_data_and_route() {
    let base_field = FiniteFieldDescriptor::new(43, NonZeroU32::new(1).expect("1 is non-zero"))
        .expect("prime field descriptor should build");
    let original = FrobeniusTrace::from_order(base_field.clone(), 52)
        .expect("original Frobenius package should build");
    let twist =
        FrobeniusTrace::from_order(base_field, 36).expect("twist Frobenius package should build");
    let point_order_report = f7_curve()
        .point_order_from_multiple(
            &f7_curve()
                .point(F7::from_i64(2), F7::from_i64(1))
                .expect("sample point should lie on the curve"),
            BigUint::from(6u8),
            &[(BigUint::from(2u8), 1), (BigUint::from(3u8), 1)],
        )
        .expect("known-multiple route should recover a sample order");
    let mestre_report = MestreGroupOrderReport::new(
        MestreConfig::with_iteration_cap(8),
        MestreSide::QuadraticTwist,
        original.clone(),
        twist.clone(),
        vec![MestreStepReport::new(
            MestreSide::QuadraticTwist,
            45,
            point_order_report,
            BigUint::from(9u8),
        )],
    );
    let report = GroupOrderReport::MestreFp(Box::new(mestre_report));

    assert_eq!(
        report.strategy(),
        GroupOrderStrategy::MestreFp(MestreConfig::with_iteration_cap(8))
    );
    assert_eq!(report.field_order(), 43);
    assert_eq!(report.curve_order(), 52);
    assert_eq!(report.trace(), -8);
    assert_eq!(report.hasse_interval(), original.hasse_interval());

    let GroupOrderReport::MestreFp(mestre_report) = report else {
        panic!("Mestre report should preserve its own variant");
    };

    assert_eq!(mestre_report.resolved_side(), MestreSide::QuadraticTwist);
    assert_eq!(mestre_report.curve_order(), 52);
    assert_eq!(mestre_report.twist_curve_order(), 36);
    assert_eq!(
        mestre_report.original_exponent_lower_bound(),
        BigUint::from(1u8)
    );
    assert_eq!(
        mestre_report.twist_exponent_lower_bound(),
        BigUint::from(9u8)
    );
    assert_eq!(mestre_report.steps().len(), 1);
    assert_eq!(mestre_report.steps()[0].side(), MestreSide::QuadraticTwist);
    assert_eq!(mestre_report.steps()[0].annihilating_multiple(), 45);
}

#[test]
fn mestre_strategy_is_honestly_reported_as_not_implemented_yet() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    assert_eq!(
        curve.group_order_by(GroupOrderStrategy::MestreFp(MestreConfig::unbounded())),
        Err(CurveError::GroupOrderStrategyRequiresSampler {
            strategy: "MestreFp"
        })
    );
}

#[test]
fn curve_type_report_classifies_an_ordinary_f43_curve() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");
    let report = trace.curve_type_report();

    assert_eq!(trace.trace(), 8);
    assert_eq!(trace.curve_type(), FrobeniusCurveType::Ordinary);
    assert_eq!(report.frobenius_trace(), &trace);
    assert_eq!(report.curve_type(), FrobeniusCurveType::Ordinary);
    assert_eq!(report.trace_mod_characteristic(), 8);
    assert!(!report.characteristic_divides_trace());
    assert!(report.is_ordinary());
    assert!(!report.is_supersingular());
}

#[test]
fn curve_type_report_classifies_a_supersingular_f43_curve() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::zero()).expect("valid curve");
    let trace = curve
        .frobenius_trace()
        .expect("Frobenius trace should compute");
    let report = trace.curve_type_report();

    assert_eq!(trace.trace(), 0);
    assert_eq!(trace.curve_type(), FrobeniusCurveType::Supersingular);
    assert_eq!(report.curve_type(), FrobeniusCurveType::Supersingular);
    assert_eq!(report.trace_mod_characteristic(), 0);
    assert!(report.characteristic_divides_trace());
    assert!(!report.is_ordinary());
    assert!(report.is_supersingular());
}

#[test]
fn relative_frobenius_on_exact_torsion_reports_fixed_f43_two_torsion() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid curve");

    let report = relative_frobenius_on_exact_torsion(&curve, 2)
        .expect("relative Frobenius on exact torsion should evaluate");

    assert_eq!(report.exact_order(), 2);
    assert_eq!(report.points().len(), 3);
    assert!(report.all_fixed());
    assert_eq!(report.fixed_count(), 3);
    assert_eq!(report.moved_count(), 0);
    assert_eq!(report.prime_field_rational_count(), 3);
    assert_eq!(report.extension_only_count(), 0);
    assert_eq!(report.orbit_count(), 3);
    assert_eq!(report.orbit_periods(), vec![1, 1, 1]);
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(1), None);
    assert_eq!(
        report.count_with_minimal_absolute_frobenius_fixing_power(1),
        None
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.fixed_by_frobenius())
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.point() == point.frobenius_image())
    );
}

#[test]
fn relative_frobenius_on_exact_torsion_rejects_zero_order() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::zero(), F43::one()).expect("valid curve");

    assert_eq!(
        relative_frobenius_on_exact_torsion(&curve, 0),
        Err(crate::elliptic_curves::CurveError::InvalidTorsionOrder { order: 0 })
    );
}

#[test]
fn absolute_frobenius_on_exact_four_torsion_can_move_extension_field_points() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::zero()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let report = absolute_frobenius_on_exact_torsion(&curve, 4, 1)
        .expect("absolute Frobenius on exact torsion should evaluate");

    assert_eq!(report.exact_order(), 4);
    assert_eq!(report.points().len(), 12);
    assert!(!report.all_fixed());
    assert_eq!(report.fixed_count(), 0);
    assert_eq!(report.moved_count(), 12);
    assert_eq!(report.prime_field_rational_count(), 0);
    assert_eq!(report.extension_only_count(), 12);
    assert_eq!(report.orbit_count(), 6);
    assert_eq!(report.orbit_periods(), vec![2, 2, 2, 2, 2, 2]);
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(1), Some(0));
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(2), Some(12));
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(3), Some(0));
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(4), Some(12));
    assert_eq!(
        report.count_with_minimal_absolute_frobenius_fixing_power(1),
        Some(0)
    );
    assert_eq!(
        report.count_with_minimal_absolute_frobenius_fixing_power(2),
        Some(12)
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| !point.fixed_by_frobenius())
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.minimal_absolute_frobenius_fixing_power() == Some(2))
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.fixed_by_absolute_frobenius_power(1) == Some(false))
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.fixed_by_absolute_frobenius_power(2) == Some(true))
    );
}

#[test]
fn absolute_frobenius_squared_fixes_exact_four_torsion_over_f43_squared() {
    let curve = ShortWeierstrassCurve::<F43Sqrt2>::new(
        F43Sqrt2::from_base(F43::zero()),
        F43Sqrt2::from_base(F43::one()),
    )
    .expect("base-defined curve should stay smooth over F43^2");

    let report = absolute_frobenius_on_exact_torsion(&curve, 4, 2)
        .expect("absolute Frobenius square should evaluate");

    assert_eq!(report.exact_order(), 4);
    assert_eq!(report.points().len(), 12);
    assert!(report.all_fixed());
    assert_eq!(report.fixed_count(), 12);
    assert_eq!(report.moved_count(), 0);
    assert_eq!(report.prime_field_rational_count(), 12);
    assert_eq!(report.extension_only_count(), 0);
    assert_eq!(report.orbit_count(), 12);
    assert_eq!(report.orbit_periods(), vec![1; 12]);
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(1), Some(0));
    assert_eq!(report.fixed_by_absolute_frobenius_power_count(2), Some(12));
    assert_eq!(
        report.count_with_minimal_absolute_frobenius_fixing_power(2),
        Some(12)
    );
    assert!(
        report
            .points()
            .iter()
            .all(|point| point.fixed_by_frobenius())
    );
}

#[test]
fn characteristic_equation_holds_on_a_prime_field_rational_point() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) should lie on the curve");
    let characteristic_polynomial = curve
        .frobenius_trace()
        .expect("trace should compute")
        .characteristic_polynomial();

    let check = verify_frobenius_characteristic_equation_at_point(
        &curve,
        &point,
        &characteristic_polynomial,
    )
    .expect("characteristic equation check should evaluate");

    assert_eq!(check.point(), &point);
    assert_eq!(check.pi_q(), &point);
    assert_eq!(check.pi_q_squared(), &point);
    assert!(check.holds());
    assert!(curve.is_identity(check.lhs()));
}

#[test]
fn characteristic_equation_rejects_a_polynomial_from_the_wrong_base_field() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) should lie on the curve");
    let wrong_polynomial = FrobeniusCharacteristicPolynomial::new(
        FiniteFieldDescriptor::new(43, core::num::NonZeroU32::new(2).expect("2 is non-zero"))
            .expect("F43^2 descriptor should be valid"),
        0,
    );

    assert_eq!(
        verify_frobenius_characteristic_equation_at_point(&curve, &point, &wrong_polynomial),
        Err(
            crate::elliptic_curves::CurveError::IncompatibleFrobeniusBaseField {
                curve_characteristic: 43,
                curve_extension_degree: 1,
                polynomial_characteristic: 43,
                polynomial_extension_degree: 2,
            }
        )
    );
}

#[test]
fn characteristic_equation_rejects_off_curve_points() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let off_curve = AffinePoint::<F43>::new(F43::zero(), F43::zero());
    let characteristic_polynomial = curve
        .frobenius_trace()
        .expect("trace should compute")
        .characteristic_polynomial();

    assert_eq!(
        verify_frobenius_characteristic_equation_at_point(
            &curve,
            &off_curve,
            &characteristic_polynomial,
        ),
        Err(crate::elliptic_curves::CurveError::PointNotOnCurve)
    );
}

#[test]
fn exhaustive_characteristic_equation_report_holds_on_a_small_prime_field_curve() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");

    let report = verify_frobenius_characteristic_equation_exhaustive(&curve)
        .expect("exhaustive characteristic-equation check should evaluate");

    assert_eq!(report.frobenius_trace().curve_order(), curve.order() as u64);
    assert_eq!(report.checked_points(), curve.order());
    assert!(report.failed_checks().is_empty());
    assert!(report.all_hold());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn property_characteristic_equation_holds_for_sampled_f43_rational_points(
        (curve, point) in arb_curve_and_point::<43>(CurveStrategyConfig::default()),
    ) {
        let characteristic_polynomial = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves")
            .characteristic_polynomial();

        let check = verify_frobenius_characteristic_equation_at_point(
            &curve,
            &point,
            &characteristic_polynomial,
        )
        .expect("characteristic equation should evaluate on enumerated rational points");

        prop_assert!(check.holds());
        prop_assert!(curve.is_identity(check.lhs()));
    }

    #[test]
    fn property_hasse_bound_holds_for_sampled_f43_curves(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
    ) {
        let report = verify_hasse_bound(&curve)
            .expect("Hasse bound should verify on small F43 curves");

        prop_assert!(report.holds());
        prop_assert!(report.trace_square() <= report.bound_square());
        prop_assert_eq!(
            report.frobenius_trace().curve_order(),
            curve.order() as u64,
        );
    }

    #[test]
    fn property_curve_type_matches_trace_divisibility_over_f43(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves");
        let report = trace.curve_type_report();
        let characteristic_divides_trace = trace.trace().rem_euclid(43) == 0;

        prop_assert_eq!(report.characteristic_divides_trace(), characteristic_divides_trace);
        prop_assert_eq!(report.trace_mod_characteristic() == 0, characteristic_divides_trace);
        prop_assert_eq!(
            report.curve_type(),
            if characteristic_divides_trace {
                FrobeniusCurveType::Supersingular
            } else {
                FrobeniusCurveType::Ordinary
            }
        );
    }

    #[test]
    fn property_relative_frobenius_on_exact_two_torsion_is_tautologically_fixed_over_f43(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
    ) {
        let report = relative_frobenius_on_exact_torsion(&curve, 2)
            .expect("relative Frobenius on exact two-torsion should evaluate");

        prop_assert!(report.all_fixed());
        prop_assert_eq!(report.fixed_count(), report.points().len());
        prop_assert_eq!(report.moved_count(), 0);
        for point in report.points() {
            prop_assert!(point.fixed_by_frobenius());
            prop_assert_eq!(point.point(), point.frobenius_image());
        }
    }

    #[test]
    fn property_exhaustive_report_matches_pointwise_checks(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
    ) {
        let report = verify_frobenius_characteristic_equation_exhaustive(&curve)
            .expect("exhaustive report should compute on small F43 curves");
        let characteristic_polynomial = report
            .frobenius_trace()
            .characteristic_polynomial();

        prop_assert_eq!(report.checked_points(), curve.order());
        prop_assert_eq!(report.all_hold(), report.failed_checks().is_empty());

        for point in curve.points() {
            let pointwise = verify_frobenius_characteristic_equation_at_point(
                &curve,
                &point,
                &characteristic_polynomial,
            )
            .expect("pointwise check should evaluate on enumerated rational points");

            let failed_in_report = report
                .failed_checks()
                .iter()
                .any(|check| check.point() == &point);

            prop_assert_eq!(pointwise.holds(), !failed_in_report);
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(4))]

    #[test]
    fn property_quadratic_twist_frobenius_relation_holds_over_f43(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
    ) {
        let package = ShortWeierstrassQuadraticTwist::new(curve, f43_quadratic_twist_factor())
            .expect("a fixed F43 non-square should define a quadratic twist package");
        let relation = package
            .frobenius_relation()
            .expect("quadratic-twist Frobenius relation should compute on small F43 curves");

        prop_assert_eq!(relation.twist_kind(), package.kind());
        prop_assert!(relation.matches_twist_kind_expectation());

        match package.kind() {
            TwistKind::Quadratic => {
                prop_assert!(relation.holds());
                prop_assert!(relation.trace_negation_holds());
                prop_assert_eq!(
                    u128::from(relation.original().curve_order())
                        + u128::from(relation.twist().curve_order()),
                    relation.expected_sum(),
                );
                prop_assert_eq!(relation.sum_orders(), relation.expected_sum());
                prop_assert_eq!(relation.twist().trace(), -relation.original().trace());
            }
            TwistKind::Trivial => {
                prop_assert!(relation.same_curve_order_holds());
                prop_assert!(relation.same_trace_holds());
                prop_assert_eq!(
                    relation.original().curve_order(),
                    relation.twist().curve_order(),
                );
                prop_assert_eq!(relation.original().trace(), relation.twist().trace());
            }
        }
    }

    #[test]
    fn property_extension_counts_of_degree_one_and_two_match_enumeration(
        curve in arb_nonsingular_curve::<17>(CurveStrategyConfig::default()),
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F17 curves");
        let lifted_curve = lift_f17_curve_to_f17_squared(&curve);

        let degree_one = compare_extension_count_with_enumeration(&curve, &trace)
            .expect("degree-one comparison should compute");
        let degree_two = compare_extension_count_with_enumeration(&lifted_curve, &trace)
            .expect("degree-two comparison should compute");

        prop_assert_eq!(degree_one.relative_extension_degree(), nz(1));
        prop_assert!(degree_one.agrees());
        prop_assert_eq!(
            degree_one.frobenius_count().curve_order(),
            &BigUint::from(curve.order() as u64),
        );

        prop_assert_eq!(degree_two.relative_extension_degree(), nz(2));
        prop_assert!(degree_two.agrees());
        prop_assert_eq!(
            degree_two.frobenius_count().curve_order(),
            &BigUint::from(lifted_curve.order() as u64),
        );
    }

    #[test]
    fn property_characteristic_polynomial_and_zeta_are_consistent_with_trace(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves");
        let polynomial = trace.characteristic_polynomial();
        let zeta_from_polynomial = polynomial.local_zeta_function();
        let zeta_from_trace = trace.local_zeta_function();

        prop_assert_eq!(polynomial.base_field(), trace.base_field());
        prop_assert_eq!(polynomial.trace(), trace.trace());
        prop_assert_eq!(polynomial.field_order(), trace.field_order());
        prop_assert_eq!(polynomial.evaluate_at_integer(1), i128::from(trace.curve_order()));
        prop_assert_eq!(&zeta_from_polynomial, &zeta_from_trace);
        prop_assert_eq!(zeta_from_polynomial.base_field(), trace.base_field());
        prop_assert_eq!(zeta_from_polynomial.trace(), trace.trace());
        prop_assert_eq!(zeta_from_polynomial.field_order(), trace.field_order());
    }

    #[test]
    fn property_scalar_isogenies_preserve_curve_order_and_trace(
        curve in arb_nonsingular_curve::<43>(CurveStrategyConfig::default()),
        scalar in 1u64..=4,
    ) {
        let trace = curve
            .frobenius_trace()
            .expect("Frobenius trace should compute on small F43 curves");
        let isogeny = ScalarMultiplicationIsogeny::new(curve, scalar)
            .expect("small scalar-multiplication isogenies should build");
        let relation = verify_isogeny_frobenius_relation(&isogeny)
            .expect("Frobenius relation should compute for scalar isogenies");

        prop_assert_eq!(relation.degree(), (scalar * scalar) as usize);
        prop_assert_eq!(relation.domain().curve_order(), trace.curve_order());
        prop_assert_eq!(relation.codomain().curve_order(), trace.curve_order());
        prop_assert_eq!(relation.domain().trace(), trace.trace());
        prop_assert_eq!(relation.codomain().trace(), trace.trace());
        prop_assert!(relation.same_curve_order());
        prop_assert!(relation.same_trace());
        prop_assert!(relation.holds());
    }

    #[test]
    fn property_isogeny_graph_frobenius_report_holds_on_small_f41_graphs(
        max_depth in 0usize..=1,
        deduplicate_by_base_field_isomorphism in any::<bool>(),
    ) {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(max_depth)
            .deduplicate_by_base_field_isomorphism(deduplicate_by_base_field_isomorphism)
            .build()
            .expect("small F41 isogeny graphs should build");
        let report = verify_isogeny_graph_frobenius_relation(&graph)
            .expect("graph Frobenius relation should compute on small F41 graphs");

        prop_assert_eq!(report.checked_nodes(), graph.node_count());
        prop_assert_eq!(report.checked_edges(), graph.edge_count());
        prop_assert!(report.all_same_curve_order());
        prop_assert!(report.all_same_trace());
        prop_assert!(report.holds());
        prop_assert!(report.node_reports().iter().all(|node| node.holds()));
    }

    #[test]
    fn property_absolute_frobenius_orbits_partition_extension_points(
        curve in arb_nonsingular_curve::<17>(CurveStrategyConfig::default()),
    ) {
        let lifted_curve = lift_f17_curve_to_f17_squared(&curve);
        let orbits = absolute_frobenius_orbits_on_points(&lifted_curve, 1)
            .expect("absolute Frobenius orbits should compute on lifted F17 curves");
        let mut seen_points = HashSet::new();
        let fixed_points_from_orbits = orbits
            .iter()
            .filter(|orbit| orbit.period() == 1)
            .count();

        for orbit in &orbits {
            prop_assert!(matches!(orbit.period(), 1 | 2));
            prop_assert_eq!(orbit.points().len(), orbit.period());
            prop_assert_eq!(
                absolute_frobenius_power_point(&lifted_curve, orbit.start(), orbit.period() as u32)
                    .expect("absolute Frobenius should evaluate on orbit starts"),
                orbit.start().clone(),
            );

            for point in orbit.points() {
                prop_assert!(seen_points.insert(point.clone()));
            }
        }

        prop_assert_eq!(
            orbits.iter().map(|orbit| orbit.period()).sum::<usize>(),
            lifted_curve.order(),
        );
        prop_assert_eq!(seen_points.len(), lifted_curve.order());
        prop_assert_eq!(fixed_points_from_orbits, curve.order());
    }
}
