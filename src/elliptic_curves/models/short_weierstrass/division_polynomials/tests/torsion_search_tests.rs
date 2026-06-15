use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve, traits::EnumerableCurveModel};
use crate::fields::{Fp, traits::Field};

use super::shared::{same_point_set, same_x_set, unique_x_coordinates_of_rational_n_torsion};

type F23 = Fp<23>;

#[test]
fn rational_roots_match_three_torsion_x_coordinates() {
    let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
        .expect("curve should be non-singular");

    let roots = curve.rational_roots_of_odd_division_polynomial(3).unwrap();
    let expected = unique_x_coordinates_of_rational_n_torsion(&curve, 3);
    assert!(same_x_set::<F23>(&roots, &expected));
}

#[test]
fn public_x_candidates_dispatch_by_parity() {
    let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
        .expect("curve should be non-singular");

    let odd = curve
        .rational_x_candidates_for_division_polynomial(3)
        .unwrap();
    let odd_expected = curve.rational_roots_of_odd_division_polynomial(3).unwrap();
    let even = curve
        .rational_x_candidates_for_division_polynomial(4)
        .unwrap();
    let even_expected = curve
        .rational_x_candidates_from_even_division_polynomial(4)
        .unwrap();

    assert!(same_x_set::<F23>(&odd, &odd_expected));
    assert!(same_x_set::<F23>(&even, &even_expected));
}

#[test]
fn public_torsion_candidates_dispatch_by_parity() {
    let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
        .expect("curve should be non-singular");

    let odd = curve
        .torsion_candidates_from_division_polynomial(3)
        .unwrap();
    let odd_expected = curve
        .torsion_candidates_from_odd_division_polynomial(3)
        .unwrap();
    let even = curve
        .torsion_candidates_from_division_polynomial(6)
        .unwrap();
    let even_expected = curve
        .torsion_candidates_from_even_division_polynomial(6)
        .unwrap();

    assert!(same_point_set::<F23>(&odd, &odd_expected));
    assert!(same_point_set::<F23>(&even, &even_expected));
}

#[test]
fn even_candidate_helpers_match_pointwise_vanishing() {
    let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
        .expect("curve should be non-singular");

    let xs = curve
        .rational_x_candidates_from_even_division_polynomial(4)
        .unwrap();
    let candidates = curve
        .torsion_candidates_from_even_division_polynomial(6)
        .unwrap();

    let mut expected_xs = Vec::new();
    for point in curve.points() {
        let AffinePoint::Finite { x, .. } = point else {
            continue;
        };
        if (curve.point_has_zero_y(&point) || curve.x_criterion_vanishes(4, &x).unwrap())
            && !expected_xs.iter().any(|existing| F23::eq(existing, &x))
        {
            expected_xs.push(x);
        }
    }

    assert!(same_x_set::<F23>(&xs, &expected_xs));
    assert!(candidates.iter().all(|point| match point {
        AffinePoint::Finite { x, .. } => {
            curve.point_has_zero_y(point) || curve.x_criterion_vanishes(6, x).unwrap()
        }
        AffinePoint::Infinity => false,
    }));
}

#[test]
fn comparison_report_matches_enumeration_on_small_example() {
    let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
        .expect("curve should be non-singular");
    let report = curve
        .compare_division_polynomial_torsion_with_enumeration(3)
        .unwrap();

    assert_eq!(report.n(), 3);
    assert!(report.missing_from_polynomial().is_empty());
    assert!(report.extra_from_polynomial().is_empty());
    assert_eq!(
        report.exact_order_polynomial_count(),
        report.exact_order_enumerated_count()
    );
}
