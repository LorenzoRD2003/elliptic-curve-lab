use proptest::prelude::*;
use std::collections::HashSet;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    traits::{CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel},
};
use crate::isogenies::{
    error::{IsogenyError, IsogenyKernelError},
    kernel::IsogenyKernel,
};

type F7 = crate::fields::Fp7;

fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
}

fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
    AffinePoint::new(F7::from_i64(x), F7::from_i64(y))
}

#[test]
fn explicit_kernel_accepts_a_valid_small_subgroup() {
    let curve = f7_curve();
    let kernel = IsogenyKernel::new(&curve, HashSet::from([curve.identity(), f7_point(6, 0)]))
        .expect("two-torsion subgroup should be valid");

    assert_eq!(kernel.order(), 2);
    assert_eq!(kernel.degree(), 2);
    assert!(kernel.contains(&curve.identity()));
    assert!(kernel.contains(&f7_point(6, 0)));
    assert_eq!(kernel.points(), &[curve.identity(), f7_point(6, 0)]);
}

#[test]
fn cyclic_kernel_enumerates_successive_multiples_of_the_generator() {
    let curve = f7_curve();
    let generator = f7_point(2, 1);
    let kernel = IsogenyKernel::cyclic(&curve, &generator)
        .expect("generator should define a cyclic subgroup");

    let expected = vec![
        curve.identity(),
        generator.clone(),
        curve
            .mul_scalar(&generator, 2)
            .expect("multiple should exist"),
        curve
            .mul_scalar(&generator, 3)
            .expect("multiple should exist"),
        curve
            .mul_scalar(&generator, 4)
            .expect("multiple should exist"),
        curve
            .mul_scalar(&generator, 5)
            .expect("multiple should exist"),
    ];

    assert_eq!(kernel.points(), expected.as_slice());
    assert_eq!(kernel.order(), 6);
    assert_eq!(kernel.degree(), 6);
}

#[test]
fn cyclic_kernel_of_identity_is_the_trivial_subgroup() {
    let curve = f7_curve();
    let kernel = IsogenyKernel::cyclic(&curve, &curve.identity())
        .expect("identity should generate the trivial subgroup");

    assert_eq!(kernel.points(), &[curve.identity()]);
    assert_eq!(kernel.order(), 1);
}

#[test]
fn constructor_rejects_empty_kernels() {
    let curve = f7_curve();

    assert!(matches!(
        IsogenyKernel::<ShortWeierstrassCurve<F7>>::new(&curve, HashSet::new()),
        Err(IsogenyError::Kernel(IsogenyKernelError::EmptyKernel))
    ));
}

#[test]
fn constructor_rejects_kernels_without_identity() {
    let curve = f7_curve();

    assert!(matches!(
        IsogenyKernel::new(&curve, HashSet::from([f7_point(6, 0)])),
        Err(IsogenyError::Kernel(
            IsogenyKernelError::KernelDoesNotContainIdentity
        ))
    ));
}

#[test]
fn constructor_rejects_points_outside_the_curve() {
    let curve = f7_curve();
    let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

    assert!(matches!(
        IsogenyKernel::new(&curve, HashSet::from([curve.identity(), invalid])),
        Err(IsogenyError::Kernel(
            IsogenyKernelError::KernelPointNotOnCurve
        ))
    ));
}

#[test]
fn constructor_rejects_sets_not_closed_under_negation() {
    let curve = f7_curve();
    let point = f7_point(2, 1);

    assert!(matches!(
        IsogenyKernel::new(&curve, HashSet::from([curve.identity(), point])),
        Err(IsogenyError::Kernel(
            IsogenyKernelError::KernelNotClosedUnderNegation
        ))
    ));
}

#[test]
fn constructor_rejects_sets_not_closed_under_addition() {
    let curve = f7_curve();
    let point = f7_point(2, 1);
    let negated = curve.neg(&point);

    assert!(matches!(
        IsogenyKernel::new(&curve, HashSet::from([curve.identity(), point, negated])),
        Err(IsogenyError::Kernel(
            IsogenyKernelError::KernelNotClosedUnderAddition
        ))
    ));
}

#[test]
fn cyclic_kernel_rejects_off_curve_generators() {
    let curve = f7_curve();
    let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

    assert!(matches!(
        IsogenyKernel::cyclic(&curve, &invalid),
        Err(IsogenyError::Curve(CurveError::PointNotOnCurve))
    ));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn property_cyclic_kernel_enumerates_exactly_the_multiples_of_the_generator(
        index in 0usize..6,
    ) {
        let curve = f7_curve();
        let points = curve.points();
        let generator = points[index].clone();
        let kernel = IsogenyKernel::cyclic(&curve, &generator)
            .expect("enumerated curve point should generate a cyclic subgroup");
        let order = curve.point_order(&generator).expect("enumerated point should have an order");
        let expected = (0..order)
            .map(|multiple| {
                curve
                    .mul_scalar(&generator, multiple as u64)
                    .expect("scalar multiples should exist")
            })
            .collect::<Vec<_>>();

        prop_assert_eq!(kernel.points(), expected.as_slice());
        prop_assert_eq!(kernel.order(), order);
    }
}
