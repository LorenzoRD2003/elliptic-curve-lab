use crate::elliptic_curves::short_weierstrass::isogenies::frobenius::{
    AbsoluteFrobeniusIsogeny, FrobeniusLikeIsogeny,
};
use crate::fields::extension_field::define_fp_quadratic_extension;
use crate::fields::traits::*;
use crate::isogenies::traits::Isogeny;

use super::shared::{F17, f17_curve};

define_fp_quadratic_extension!(
    spec: F17Sqrt3FunctionFieldSpec,
    field: F17Sqrt3Function,
        base: F17,
    non_residue: 3,
    name: "F17(sqrt(3)) for function-field Frobenius tests",
);

fn f17_sqrt3_curve() -> crate::elliptic_curves::ShortWeierstrassCurve<F17Sqrt3Function> {
    let alpha = F17Sqrt3Function::element(vec![F17::zero(), F17::one()]);
    crate::elliptic_curves::ShortWeierstrassCurve::<F17Sqrt3Function>::new(
        alpha,
        F17Sqrt3Function::one(),
    )
    .expect("curve should be nonsingular")
}

#[test]
fn inverse_absolute_frobenius_pullback_recovers_x_generator_on_the_twist() {
    let curve = f17_curve();
    let frobenius =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let twist_field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(frobenius.codomain().clone());

    let recovered = frobenius
        .x_pullback()
        .inverse_absolute_frobenius_pullback_to_twist(frobenius.codomain())
        .expect("x^p should come from the twist x-coordinate");

    assert_eq!(recovered, twist_field.x());
}

#[test]
fn inverse_absolute_frobenius_pullback_recovers_y_generator_on_the_twist() {
    let curve = f17_curve();
    let frobenius =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let twist_field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(frobenius.codomain().clone());

    let recovered = frobenius
        .y_pullback()
        .inverse_absolute_frobenius_pullback_to_twist(frobenius.codomain())
        .expect("y^p should come from the twist y-coordinate");

    assert_eq!(recovered, twist_field.y());
}

#[test]
fn inverse_absolute_frobenius_pullback_recovers_mixed_example_over_nontrivial_twist() {
    let curve = f17_sqrt3_curve();
    let frobenius =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let twist_field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17Sqrt3Function,
        >::new(frobenius.codomain().clone());
    let preimage = twist_field
        .x()
        .add(&twist_field.y())
        .expect("same-curve addition should work");
    let image = frobenius
        .as_function_field_map()
        .pullback_function(&preimage)
        .expect("Frobenius pullback should work");

    let recovered = image
        .inverse_absolute_frobenius_pullback_to_twist(frobenius.codomain())
        .expect("mixed Frobenius image should invert");

    assert_eq!(recovered, preimage);
    assert_eq!(recovered.curve(), frobenius.codomain());
    assert_ne!(frobenius.codomain(), &curve);
}

#[test]
fn inverse_absolute_frobenius_pullback_rejects_y_without_rhs_factor() {
    let curve = f17_curve();
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(curve.clone());
    let twist = curve
        .frobenius_twist_power(1)
        .expect("Frobenius twist should build");

    assert_eq!(
        field
            .y()
            .inverse_absolute_frobenius_pullback_to_twist(&twist),
        None
    );
}
