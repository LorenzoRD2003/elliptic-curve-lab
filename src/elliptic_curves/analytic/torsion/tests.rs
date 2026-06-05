use num_complex::Complex64;
use proptest::prelude::*;

use super::{
    AnalyticDivisionPolynomialComparisonCase, AnalyticDivisionPolynomialComparisonStatus,
    TorusTorsionIndex, compare_analytic_torsion_with_division_polynomial,
    compare_primitive_analytic_torsion_with_division_polynomial,
    map_primitive_torus_torsion_to_curve, map_torus_torsion_to_curve,
    primitive_torus_n_torsion_points, torus_n_torsion_points,
};
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticCurvePoint, ApproxTolerance, ComplexLattice,
    EllipticFunctionTruncation, EvenDivisionPolynomialVanishingBranch,
    FundamentalParallelogramCoordinate, LatticeSumTruncation, UpperHalfPlanePoint,
    map_torus_point_to_curve,
};

fn square_lattice() -> ComplexLattice {
    ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
}

#[test]
fn torus_torsion_index_validates_reduced_bounds() {
    let index = TorusTorsionIndex::new(1, 2, 3).expect("reduced index should be valid");

    assert_eq!(index.a(), 1);
    assert_eq!(index.b(), 2);
    assert_eq!(index.n(), 3);
    assert_eq!(
        TorusTorsionIndex::new(0, 0, 0),
        Err(AnalyticCurveError::InvalidTorusTorsionIndex)
    );
    assert_eq!(
        TorusTorsionIndex::new(3, 0, 3),
        Err(AnalyticCurveError::InvalidTorusTorsionIndex)
    );
    assert_eq!(
        TorusTorsionIndex::new(0, 3, 3),
        Err(AnalyticCurveError::InvalidTorusTorsionIndex)
    );
}

#[test]
fn primitive_index_uses_gcd_of_a_b_and_n() {
    assert!(
        TorusTorsionIndex::new(1, 1, 2)
            .expect("reduced index should be valid")
            .is_primitive()
    );
    assert!(
        !TorusTorsionIndex::new(2, 0, 4)
            .expect("reduced index should be valid")
            .is_primitive()
    );
    assert!(
        TorusTorsionIndex::new(0, 0, 1)
            .expect("identity in one torsion should be primitive by convention")
            .is_primitive()
    );
}

#[test]
fn identity_helpers_detect_only_the_zero_zero_class() {
    let identity = TorusTorsionIndex::new(0, 0, 3).unwrap();
    let non_identity = TorusTorsionIndex::new(0, 1, 3).unwrap();

    assert!(identity.is_identity_class());
    assert!(!non_identity.is_identity_class());

    let lattice = square_lattice();
    let points = torus_n_torsion_points(&lattice, 2).unwrap();
    assert!(points[0].is_identity_class());
    assert!(!points[1].is_identity_class());
}

#[test]
fn torus_n_torsion_points_reject_zero_order() {
    let lattice = square_lattice();

    assert_eq!(
        torus_n_torsion_points(&lattice, 0),
        Err(AnalyticCurveError::InvalidTorusTorsionIndex)
    );
}

#[test]
fn torus_n_torsion_points_return_n_squared_points_in_lexicographic_order() {
    let lattice = square_lattice();
    let points = torus_n_torsion_points(&lattice, 2).expect("two torsion should succeed");

    assert_eq!(points.len(), 4);
    assert_eq!(points[0].index().a(), 0);
    assert_eq!(points[0].index().b(), 0);
    assert_eq!(points[1].index().a(), 0);
    assert_eq!(points[1].index().b(), 1);
    assert_eq!(points[2].index().a(), 1);
    assert_eq!(points[2].index().b(), 0);
    assert_eq!(points[3].index().a(), 1);
    assert_eq!(points[3].index().b(), 1);
}

#[test]
fn two_torsion_points_match_expected_square_lattice_representatives() {
    let lattice = square_lattice();
    let points = torus_n_torsion_points(&lattice, 2).expect("two torsion should succeed");

    assert_eq!(
        points[0].coordinate(),
        &FundamentalParallelogramCoordinate::new(0.0, 0.0).unwrap()
    );
    assert_eq!(
        points[1].coordinate(),
        &FundamentalParallelogramCoordinate::new(0.0, 0.5).unwrap()
    );
    assert_eq!(
        points[2].coordinate(),
        &FundamentalParallelogramCoordinate::new(0.5, 0.0).unwrap()
    );
    assert_eq!(
        points[3].coordinate(),
        &FundamentalParallelogramCoordinate::new(0.5, 0.5).unwrap()
    );
    assert_eq!(points[0].z(), &Complex64::new(0.0, 0.0));
    assert_eq!(points[1].z(), &Complex64::new(0.0, 0.5));
    assert_eq!(points[2].z(), &Complex64::new(0.5, 0.0));
    assert_eq!(points[3].z(), &Complex64::new(0.5, 0.5));
}

#[test]
fn primitive_two_torsion_excludes_only_the_identity() {
    let lattice = square_lattice();
    let points = primitive_torus_n_torsion_points(&lattice, 2)
        .expect("primitive two torsion should succeed");

    assert_eq!(points.len(), 3);
    assert!(points.iter().all(|point| point.index().is_primitive()));
    assert!(
        !points
            .iter()
            .any(|point| point.index().a() == 0 && point.index().b() == 0)
    );
}

#[test]
fn primitive_three_torsion_has_eight_points() {
    let lattice = square_lattice();
    let points = primitive_torus_n_torsion_points(&lattice, 3)
        .expect("primitive three torsion should succeed");

    assert_eq!(points.len(), 8);
    assert!(points.iter().all(|point| point.index().is_primitive()));
}

#[test]
fn one_torsion_returns_the_identity_class() {
    let lattice = square_lattice();
    let points = torus_n_torsion_points(&lattice, 1).expect("one torsion should succeed");
    let primitive = primitive_torus_n_torsion_points(&lattice, 1)
        .expect("primitive one torsion should succeed");

    assert_eq!(points.len(), 1);
    assert_eq!(primitive.len(), 1);
    assert_eq!(points[0].z(), &Complex64::new(0.0, 0.0));
    assert_eq!(primitive[0].index(), points[0].index());
}

#[test]
fn mapping_torus_torsion_rejects_zero_order() {
    let lattice = square_lattice();

    assert_eq!(
        map_torus_torsion_to_curve(
            &lattice,
            0,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::strict(),
        ),
        Err(AnalyticCurveError::InvalidTorusTorsionIndex)
    );
}

#[test]
fn identity_torsion_class_maps_to_infinity() {
    let lattice = square_lattice();
    let mapped = map_torus_torsion_to_curve(
        &lattice,
        2,
        LatticeSumTruncation::default_educational(),
        EllipticFunctionTruncation::default_educational(),
        ApproxTolerance::strict(),
    )
    .expect("two torsion mapping should succeed");

    let identity = &mapped[0];
    assert_eq!(identity.torus_point().index().a(), 0);
    assert_eq!(identity.torus_point().index().b(), 0);
    assert_eq!(identity.curve_point(), &AnalyticCurvePoint::infinity());
    assert!(identity.lies_on_curve());
}

#[test]
fn mapped_torsion_points_preserve_lexicographic_order_and_count() {
    let lattice = square_lattice();
    let mapped = map_torus_torsion_to_curve(
        &lattice,
        2,
        LatticeSumTruncation::default_educational(),
        EllipticFunctionTruncation::default_educational(),
        ApproxTolerance::loose(),
    )
    .expect("two torsion mapping should succeed");

    assert_eq!(mapped.len(), 4);
    assert_eq!(mapped[0].torus_point().index().a(), 0);
    assert_eq!(mapped[0].torus_point().index().b(), 0);
    assert_eq!(mapped[1].torus_point().index().a(), 0);
    assert_eq!(mapped[1].torus_point().index().b(), 1);
    assert_eq!(mapped[2].torus_point().index().a(), 1);
    assert_eq!(mapped[2].torus_point().index().b(), 0);
    assert_eq!(mapped[3].torus_point().index().a(), 1);
    assert_eq!(mapped[3].torus_point().index().b(), 1);
}

#[test]
fn mapped_torsion_points_match_pointwise_torus_to_curve_mapping() {
    let lattice = square_lattice();
    let invariant_truncation = LatticeSumTruncation::larger_for_comparison();
    let function_truncation = EllipticFunctionTruncation::default_educational();
    let tolerance = ApproxTolerance::loose();
    let mapped = map_torus_torsion_to_curve(
        &lattice,
        2,
        invariant_truncation,
        function_truncation,
        tolerance,
    )
    .expect("two torsion mapping should succeed");

    for entry in mapped {
        let direct = map_torus_point_to_curve(
            &lattice,
            *entry.torus_point().z(),
            invariant_truncation,
            function_truncation,
            tolerance,
        )
        .expect("pointwise torus-to-curve map should succeed");

        assert_eq!(entry.curve_point(), direct.point());
        assert_eq!(entry.membership_report(), direct.membership_report());
    }
}

#[test]
fn mapped_primitive_two_torsion_excludes_only_the_identity() {
    let lattice = square_lattice();
    let mapped = map_primitive_torus_torsion_to_curve(
        &lattice,
        2,
        LatticeSumTruncation::default_educational(),
        EllipticFunctionTruncation::default_educational(),
        ApproxTolerance::loose(),
    )
    .expect("primitive two torsion mapping should succeed");

    assert_eq!(mapped.len(), 3);
    assert!(
        !mapped
            .iter()
            .any(|point| point.torus_point().index().a() == 0
                && point.torus_point().index().b() == 0)
    );
    assert!(
        mapped
            .iter()
            .all(|point| point.torus_point().index().is_primitive())
    );
}

#[test]
fn mapped_primitive_one_torsion_keeps_the_identity_infinity_point() {
    let lattice = square_lattice();
    let mapped = map_primitive_torus_torsion_to_curve(
        &lattice,
        1,
        LatticeSumTruncation::default_educational(),
        EllipticFunctionTruncation::default_educational(),
        ApproxTolerance::strict(),
    )
    .expect("primitive one torsion mapping should succeed");

    assert_eq!(mapped.len(), 1);
    assert_eq!(mapped[0].curve_point(), &AnalyticCurvePoint::infinity());
    assert!(mapped[0].torus_point().index().is_primitive());
}

#[test]
fn analytic_division_polynomial_comparison_reports_pole_at_identity() {
    let lattice = square_lattice();
    let reports = compare_analytic_torsion_with_division_polynomial(
        &lattice,
        2,
        LatticeSumTruncation::default_educational(),
        EllipticFunctionTruncation::default_educational(),
        ApproxTolerance::loose(),
    )
    .expect("comparison should succeed");

    let identity = &reports[0];
    match identity {
        AnalyticDivisionPolynomialComparisonCase::Pole {
            torsion_point,
            tolerance,
        } => {
            assert_eq!(torsion_point.curve_point(), &AnalyticCurvePoint::infinity());
            assert_eq!(*tolerance, ApproxTolerance::loose());
        }
        other => panic!("expected pole report, got {other:?}"),
    }
}

#[test]
fn analytic_division_polynomial_comparison_reports_even_branch_state_for_two_torsion() {
    let lattice = square_lattice();
    let reports = compare_primitive_analytic_torsion_with_division_polynomial(
        &lattice,
        2,
        LatticeSumTruncation::default_educational(),
        EllipticFunctionTruncation::default_educational(),
        ApproxTolerance::loose(),
    )
    .expect("primitive two-torsion comparison should succeed");

    assert_eq!(reports.len(), 3);
    assert!(
        reports
            .iter()
            .all(|report| matches!(report, AnalyticDivisionPolynomialComparisonCase::Even(_)))
    );
    assert!(reports.iter().all(|report| match report {
        AnalyticDivisionPolynomialComparisonCase::Even(even_report) => {
            even_report.branch() == &EvenDivisionPolynomialVanishingBranch::NeitherBranch
                && even_report.status()
                    == &AnalyticDivisionPolynomialComparisonStatus::DoesNotVanishApproximately
        }
        _ => false,
    }));
}

#[test]
fn analytic_division_polynomial_comparison_improves_with_larger_truncations_for_three_torsion() {
    let lattice = square_lattice();
    let small_reports = compare_primitive_analytic_torsion_with_division_polynomial(
        &lattice,
        3,
        LatticeSumTruncation::larger_for_comparison(),
        EllipticFunctionTruncation::new(6).unwrap(),
        ApproxTolerance::new(1.0e-2, 1.0e-2),
    )
    .expect("primitive three-torsion comparison should succeed");
    let large_reports = compare_primitive_analytic_torsion_with_division_polynomial(
        &lattice,
        3,
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        ApproxTolerance::new(1.0e-2, 1.0e-2),
    )
    .expect("primitive three-torsion comparison should succeed");

    assert_eq!(small_reports.len(), 8);
    assert_eq!(large_reports.len(), 8);
    assert!(
        large_reports
            .iter()
            .all(|report| matches!(report, AnalyticDivisionPolynomialComparisonCase::Odd(_)))
    );

    let small_max = small_reports
        .iter()
        .filter_map(|report| match report {
            AnalyticDivisionPolynomialComparisonCase::Odd(odd_report) => {
                Some(odd_report.absolute_value())
            }
            _ => None,
        })
        .fold(0.0_f64, f64::max);
    let large_max = large_reports
        .iter()
        .filter_map(|report| match report {
            AnalyticDivisionPolynomialComparisonCase::Odd(odd_report) => {
                Some(odd_report.absolute_value())
            }
            _ => None,
        })
        .fold(0.0_f64, f64::max);

    assert!(large_max < small_max);
}

#[test]
fn analytic_division_polynomial_comparison_rejects_zero_order() {
    let lattice = square_lattice();

    assert_eq!(
        compare_analytic_torsion_with_division_polynomial(
            &lattice,
            0,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::strict(),
        ),
        Err(AnalyticCurveError::InvalidTorusTorsionIndex)
    );
}

fn gcd_usize(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }

    a
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn valid_torus_torsion_indices_match_identity_and_primitivity_criteria(
        n in 1usize..8,
        a in 0usize..8,
        b in 0usize..8,
    ) {
        if a < n && b < n {
            let index = TorusTorsionIndex::new(a, b, n).unwrap();
            let expected_primitive = gcd_usize(gcd_usize(a, b), n) == 1;

            prop_assert_eq!(index.a(), a);
            prop_assert_eq!(index.b(), b);
            prop_assert_eq!(index.n(), n);
            prop_assert_eq!(index.is_identity_class(), a == 0 && b == 0);
            prop_assert_eq!(index.is_primitive(), expected_primitive);
        } else {
            prop_assert_eq!(
                TorusTorsionIndex::new(a, b, n),
                Err(AnalyticCurveError::InvalidTorusTorsionIndex)
            );
        }
    }

    #[test]
    fn torus_n_torsion_counts_match_the_expected_arithmetic(
        n in 1usize..8,
    ) {
        let lattice = square_lattice();
        let points = torus_n_torsion_points(&lattice, n).unwrap();
        let primitive = primitive_torus_n_torsion_points(&lattice, n).unwrap();
        let expected_primitive_count = (0..n)
            .flat_map(|a| (0..n).map(move |b| (a, b)))
            .filter(|(a, b)| gcd_usize(gcd_usize(*a, *b), n) == 1)
            .count();

        prop_assert_eq!(points.len(), n * n);
        prop_assert_eq!(primitive.len(), expected_primitive_count);
        prop_assert!(points[0].is_identity_class());
        prop_assert!(primitive.iter().all(|point| point.index().is_primitive()));
    }
}
