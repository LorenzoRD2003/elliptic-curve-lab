use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{AffineCurveModel, CurveModel, GroupCurveModel};
use crate::fields::RationalFunction;
use crate::fields::{Field, Fp};
use crate::isogenies::velu::VeluIsogeny;
use crate::isogenies::{
    Isogeny, IsogenyError, IsogenyKernel, IsogenyKernelError, IsogenySeparabilityKind,
};
use crate::polynomials::evaluation::evaluate_dense;
use crate::proptest_support::isogenies::arb_cyclic_kernel_case;
use proptest::prelude::*;
use std::collections::HashSet;

type F41 = Fp<41>;

fn f41_curve() -> ShortWeierstrassCurve<F41> {
    ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

#[test]
fn from_points_rejects_invalid_kernel_before_reaching_velu_todo() {
    let domain = f41_curve();
    let invalid = domain
        .point(F41::from_i64(3), F41::from_i64(6))
        .expect("point should lie on the curve");

    let result = VeluIsogeny::from_points(domain, HashSet::from([invalid]));

    assert!(matches!(
        result,
        Err(IsogenyError::Kernel(
            IsogenyKernelError::KernelDoesNotContainIdentity
        ))
    ));
}

#[test]
fn from_generator_rejects_off_curve_points_before_reaching_velu_todo() {
    let domain = f41_curve();
    let invalid =
        crate::elliptic_curves::affine::AffinePoint::<F41>::new(F41::from_i64(2), F41::from_i64(2));

    let result = VeluIsogeny::from_generator(domain, invalid);

    assert!(matches!(result, Err(IsogenyError::Curve(_))));
}

#[test]
fn from_generator_builds_the_current_velu_scaffold_after_kernel_validation() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");

    let isogeny =
        VeluIsogeny::from_generator(domain.clone(), generator.clone()).expect("should build");

    assert_eq!(isogeny.domain().to_string(), domain.to_string());
    assert_eq!(isogeny.degree(), 2);
    assert_eq!(isogeny.kernel_points().len(), 2);
    assert_eq!(isogeny.kernel_nonzero_points(), &[generator]);
    assert_eq!(
        isogeny.codomain().to_string(),
        "y^2 = x^3 + (18 (mod 41))x + (38 (mod 41))"
    );
}

#[test]
fn evaluate_maps_the_f41_example_point_into_the_codomain() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let isogeny =
        VeluIsogeny::from_generator(domain.clone(), generator.clone()).expect("should build");
    let point = domain
        .point(F41::from_i64(3), F41::from_i64(6))
        .expect("point should lie on the curve");
    let expected = isogeny
        .codomain()
        .point(F41::from_i64(35), F41::from_i64(40))
        .expect("the translated point should lie on the codomain");

    assert_eq!(isogeny.evaluate(&point), Ok(expected));
}

#[test]
fn evaluate_sends_kernel_points_to_the_codomain_identity() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let isogeny = VeluIsogeny::from_generator(domain, generator.clone()).expect("should build");

    assert_eq!(
        isogeny
            .evaluate(&AffinePoint::infinity())
            .expect("evaluation should succeed"),
        isogeny.codomain().identity()
    );
    assert_eq!(
        isogeny
            .evaluate(&generator)
            .expect("evaluation should succeed"),
        isogeny.codomain().identity()
    );
}

#[test]
fn translation_sum_coordinates_match_the_f41_two_torsion_example() {
    let domain = f41_curve();
    let codomain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let kernel =
        IsogenyKernel::cyclic(&domain, &generator).expect("two-torsion subgroup should work");
    let isogeny = VeluIsogeny {
        degree: kernel.degree(),
        domain: domain.clone(),
        codomain,
        kernel,
    };
    let point = domain
        .point(F41::from_i64(3), F41::from_i64(6))
        .expect("point should lie on the curve");

    assert_eq!(
        isogeny
            .translation_sum_coordinates(&point)
            .expect("translation sums should succeed"),
        Some((F41::from_i64(35), F41::from_i64(40)))
    );
}

#[test]
fn velu_codomain_curve_matches_the_f41_two_torsion_example() {
    assert_eq!(
        VeluIsogeny::from_generator(
            f41_curve(),
            f41_curve()
                .point(F41::from_i64(40), F41::from_i64(0))
                .expect("point should lie on the curve")
        )
        .expect("isogeny should build")
        .codomain()
        .to_string(),
        "y^2 = x^3 + (18 (mod 41))x + (38 (mod 41))"
    );
}

#[test]
fn function_field_map_reuses_the_exported_pullback_generators() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let isogeny = VeluIsogeny::from_generator(domain, generator).expect("should build");
    let map = isogeny.as_function_field_map();

    assert_eq!(map.x_pullback(), &isogeny.x_pullback());
    assert_eq!(map.y_pullback(), &isogeny.y_pullback());
}

#[test]
fn function_field_map_domain_and_codomain_match_the_velu_isogeny() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let isogeny = VeluIsogeny::from_generator(domain.clone(), generator).expect("should build");
    let map = isogeny.as_function_field_map();

    assert_eq!(map.domain_curve(), isogeny.domain());
    assert_eq!(map.codomain_curve(), isogeny.codomain());
}

#[test]
fn function_field_pullbacks_recover_point_evaluation_away_from_the_kernel() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let isogeny =
        VeluIsogeny::from_generator(domain.clone(), generator.clone()).expect("should build");
    let point = domain
        .point(F41::from_i64(3), F41::from_i64(6))
        .expect("point should lie on the curve");
    let image = isogeny.evaluate(&point).expect("point should evaluate");
    let x_pullback = isogeny.x_pullback();
    let y_pullback = isogeny.y_pullback();
    let x_value = evaluate_short_weierstrass_function_at_point(&x_pullback, &point)
        .expect("x pullback should be regular away from the kernel");
    let y_value = evaluate_short_weierstrass_function_at_point(&y_pullback, &point)
        .expect("y pullback should be regular away from the kernel");

    assert_eq!(AffinePoint::x_coordinate(&image), Some(&x_value));
    assert_eq!(AffinePoint::y_coordinate(&image), Some(&y_value));
}

#[test]
fn velu_differential_pullback_report_matches_the_exported_function_field_map() {
    let domain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let isogeny = VeluIsogeny::from_generator(domain, generator).expect("should build");
    let direct = isogeny
        .as_function_field_map()
        .differential_pullback_report()
        .expect("map report should build");
    let wrapped = isogeny
        .differential_pullback_report()
        .expect("Velu report should build");

    assert_eq!(wrapped.rational_multiplier(), direct.rational_multiplier());
    assert_eq!(
        wrapped.separability_kind(),
        IsogenySeparabilityKind::Separable
    );
}

#[test]
fn translation_sum_coordinates_return_none_on_kernel_points() {
    let domain = f41_curve();
    let codomain = f41_curve();
    let generator = domain
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("point should lie on the curve");
    let kernel =
        IsogenyKernel::cyclic(&domain, &generator).expect("two-torsion subgroup should work");
    let isogeny = VeluIsogeny {
        degree: kernel.degree(),
        domain,
        codomain,
        kernel,
    };

    assert!(
        isogeny
            .translation_sum_coordinates(&AffinePoint::infinity())
            .expect("kernel identity should produce no affine coordinates")
            .is_none(),
    );
    assert!(
        isogeny
            .translation_sum_coordinates(&generator)
            .expect("kernel point should produce no affine coordinates")
            .is_none()
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn property_velu_is_constant_on_kernel_cosets(
        case in arb_cyclic_kernel_case(),
    ) {
        let sample_image = case
            .isogeny
            .evaluate(&case.sample_point)
            .expect("sample point should evaluate");
        let coset_image = case
            .isogeny
            .evaluate(&case.coset_point)
            .expect("kernel translate should evaluate");

        prop_assert_eq!(case.isogeny.degree(), case.isogeny.kernel_points().len());
        prop_assert_eq!(
            case.curve
                .add(&case.sample_point, &case.kernel_point)
                .expect("kernel translation should stay on the curve"),
            case.coset_point
        );
        prop_assert_eq!(sample_image, coset_image);
    }
}

fn evaluate_short_weierstrass_function_at_point<F: Field>(
    function: &crate::elliptic_curves::ShortWeierstrassFunction<F>,
    point: &AffinePoint<F>,
) -> Option<F::Elem> {
    let x = AffinePoint::x_coordinate(point)?;
    let y = AffinePoint::y_coordinate(point)?;
    let a_value = evaluate_rational_function_at_x(function.a_part(), x)?;
    let b_value = evaluate_rational_function_at_x(function.b_part(), x)?;

    Some(F::add(&a_value, &F::mul(y, &b_value)))
}

fn evaluate_rational_function_at_x<F: Field>(
    function: &RationalFunction<F>,
    x: &F::Elem,
) -> Option<F::Elem> {
    let numerator = evaluate_dense(function.numerator(), x).ok()?;
    let denominator = evaluate_dense(function.denominator(), x).ok()?;

    if F::is_zero(&denominator) {
        None
    } else {
        F::div(&numerator, &denominator).ok()
    }
}
