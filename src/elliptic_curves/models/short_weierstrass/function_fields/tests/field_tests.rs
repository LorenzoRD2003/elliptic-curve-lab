use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::traits::AmbientField;

use super::shared::{F17, f17_curve};

#[test]
fn function_field_family_exposes_zero_one_x_and_y() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());

    assert!(field.zero().is_zero());
    assert!(field.one().is_one());
    assert_eq!(
        field.x().a_part(),
        &crate::fields::rational_function_field::RationalFunction::<F17>::indeterminate()
    );
    assert!(field.x().b_part().is_zero());
    assert!(field.y().a_part().is_zero());
    assert!(field.y().b_part().is_one());
}

#[test]
fn ambient_field_zero_one_and_equality_match_the_function_field_family() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
    let zero = AmbientField::zero(&field);
    let one = AmbientField::one(&field);

    assert!(AmbientField::is_zero(&field, &zero));
    assert!(AmbientField::is_one(&field, &one));
    assert!(AmbientField::eq(&field, &field.x(), &field.x()));
}

#[test]
fn ambient_field_default_sub_and_div_work_for_function_fields() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
    let x = field.x();
    let one = field.one();
    let x_plus_one = AmbientField::add(&field, &x, &one).expect("addition should work");
    let recovered_one =
        AmbientField::sub(&field, &x_plus_one, &x).expect("default subtraction should work");
    let recovered_x = AmbientField::div(
        &field,
        &x.mul(&one).expect("multiplication should work"),
        &one,
    )
    .expect("default division should work");

    assert!(AmbientField::eq(&field, &recovered_one, &one));
    assert!(AmbientField::eq(&field, &recovered_x, &x));
}

#[test]
fn ambient_field_reports_incompatible_curve_operations() {
    let first_curve = f17_curve();
    let second_curve = ShortWeierstrassCurve::<crate::fields::Fp17>::new(
        crate::fields::Fp17::from_i64(5),
        crate::fields::Fp17::from_i64(7),
    )
    .expect("curve should be nonsingular");
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(first_curve.clone());
    let left = ShortWeierstrassFunction::<F17>::one(first_curve);
    let right = ShortWeierstrassFunction::<F17>::one(second_curve);

    assert_eq!(
        AmbientField::add(&field, &left, &right),
        Err(CurveError::IncompatibleFunctionFieldCurves)
    );
    assert_eq!(
        AmbientField::mul(&field, &left, &right),
        Err(CurveError::IncompatibleFunctionFieldCurves)
    );
}
