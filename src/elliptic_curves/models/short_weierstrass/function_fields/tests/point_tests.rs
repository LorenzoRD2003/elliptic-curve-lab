use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::function_fields::{
        ShortWeierstrassFunctionField, ShortWeierstrassFunctionFieldPoint,
    },
    traits::{AffineCurveModel, GroupCurveModel},
};
use crate::fields::{Fp, traits::Field};

use super::shared::{F17, evaluate_short_weierstrass_function_at_point, f17_curve};

#[test]
fn generic_point_is_the_pair_of_distinguished_coordinate_functions() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let generic = field.generic_point();

    assert_eq!(generic.x(), Some(&field.x()));
    assert_eq!(generic.y(), Some(&field.y()));
}

#[test]
fn translating_the_generic_point_by_infinity_is_the_identity() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let translated = field
        .translate_generic_point_by_base_point(&AffinePoint::Infinity)
        .expect("identity translation should succeed");

    assert_eq!(translated, field.generic_point());
}

#[test]
fn adding_infinity_is_neutral_for_function_field_points() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let generic = field.generic_point();

    assert_eq!(
        field.add_points(&generic, &ShortWeierstrassFunctionFieldPoint::Infinity),
        Ok(generic.clone())
    );
    assert_eq!(
        field.add_points(&ShortWeierstrassFunctionFieldPoint::Infinity, &generic),
        Ok(generic)
    );
}

#[test]
fn adding_a_point_and_its_negative_gives_infinity() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let point = field.generic_point();
    let negated = field.neg_point(&point).expect("negation should succeed");

    assert_eq!(
        field.add_points(&point, &negated),
        Ok(ShortWeierstrassFunctionFieldPoint::Infinity)
    );
}

#[test]
fn adding_a_constant_point_to_the_generic_point_matches_translation_helper() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let translation_point = curve
        .point(Fp::<41>::from_i64(40), Fp::<41>::from_i64(0))
        .expect("sample translation point should lie on the curve");
    let generic = field.generic_point();
    let constant = field
        .embed_constant_affine_point(&translation_point)
        .expect("constant embedding should succeed");

    assert_eq!(
        field.add_points(&generic, &constant),
        field.translate_generic_point_by_base_point(&translation_point)
    );
}

#[test]
fn doubling_the_generic_point_matches_the_generic_doubling_helper() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let generic = field.generic_point();

    assert_eq!(field.double_point(&generic), field.double_generic_point());
}

#[test]
fn generic_scalar_multiplication_by_three_matches_curve_scalar_multiplication_on_a_sample_point() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let sample_point = curve
        .point(Fp::<41>::from_i64(3), Fp::<41>::from_i64(6))
        .expect("sample point should lie on the curve");
    let expected = curve
        .mul_scalar(&sample_point, 3)
        .expect("scalar multiplication should stay on the curve");
    let triple = field
        .generic_point_multiple(3)
        .expect("generic scalar multiplication should succeed");
    let triple_x = triple.x().expect("triple stays affine");
    let triple_y = triple.y().expect("triple stays affine");

    assert_eq!(
        evaluate_short_weierstrass_function_at_point(triple_x, &sample_point),
        AffinePoint::x_coordinate(&expected).cloned()
    );
    assert_eq!(
        evaluate_short_weierstrass_function_at_point(triple_y, &sample_point),
        AffinePoint::y_coordinate(&expected).cloned()
    );
}

#[test]
fn translating_the_generic_point_matches_curve_addition_on_a_sample_point() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let sample_point = curve
        .point(Fp::<41>::from_i64(3), Fp::<41>::from_i64(6))
        .expect("sample point should lie on the curve");
    let translation_point = curve
        .point(Fp::<41>::from_i64(40), Fp::<41>::from_i64(0))
        .expect("sample translation point should lie on the curve");
    let expected = curve
        .add(&sample_point, &translation_point)
        .expect("sample translation should stay on the curve");
    let translated = field
        .translate_generic_point_by_base_point(&translation_point)
        .expect("generic translation should succeed");
    let translated_x = translated.x().expect("translation stays affine");
    let translated_y = translated.y().expect("translation stays affine");

    assert_eq!(
        evaluate_short_weierstrass_function_at_point(translated_x, &sample_point),
        AffinePoint::x_coordinate(&expected).cloned()
    );
    assert_eq!(
        evaluate_short_weierstrass_function_at_point(translated_y, &sample_point),
        AffinePoint::y_coordinate(&expected).cloned()
    );
}

#[test]
fn doubling_the_generic_point_matches_curve_doubling_on_a_sample_point() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let sample_point = curve
        .point(Fp::<41>::from_i64(3), Fp::<41>::from_i64(6))
        .expect("sample point should lie on the curve");
    let expected = curve
        .double(&sample_point)
        .expect("sample doubling should stay on the curve");
    let doubled = field
        .double_generic_point()
        .expect("generic doubling should succeed");
    let doubled_x = doubled.x().expect("doubling stays affine");
    let doubled_y = doubled.y().expect("doubling stays affine");

    assert_eq!(
        evaluate_short_weierstrass_function_at_point(doubled_x, &sample_point),
        AffinePoint::x_coordinate(&expected).cloned()
    );
    assert_eq!(
        evaluate_short_weierstrass_function_at_point(doubled_y, &sample_point),
        AffinePoint::y_coordinate(&expected).cloned()
    );
}
