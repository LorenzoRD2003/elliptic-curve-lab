use super::VeluIsogeny;
use crate::elliptic_curves::{AffineCurveModel, AffinePoint, CurveModel, ShortWeierstrassCurve};
use crate::fields::{Field, Fp};
use crate::isogenies::{Isogeny, IsogenyError, IsogenyKernel};
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
        Err(IsogenyError::KernelDoesNotContainIdentity)
    ));
}

#[test]
fn from_generator_rejects_off_curve_points_before_reaching_velu_todo() {
    let domain = f41_curve();
    let invalid =
        crate::elliptic_curves::AffinePoint::<F41>::new(F41::from_i64(2), F41::from_i64(2));

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
