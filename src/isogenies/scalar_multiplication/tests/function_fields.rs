use super::shared::{F41, curve, evaluate_short_weierstrass_function_at_point};
use crate::elliptic_curves::AffinePoint;
use crate::elliptic_curves::{
    short_weierstrass::{
        function_fields::ShortWeierstrassFunctionField,
        isogenies::function_field_maps::ShortWeierstrassFunctionFieldMap,
    },
    traits::EnumerableCurveModel,
};
use crate::isogenies::{scalar_multiplication::ScalarMultiplicationIsogeny, traits::Isogeny};

#[test]
fn function_field_map_of_scalar_one_is_the_identity_pullback() {
    let curve = curve();
    let field = ShortWeierstrassFunctionField::<F41>::new(curve.clone());
    let scalar = ScalarMultiplicationIsogeny::new(curve.clone(), 1)
        .expect("scalar multiplication should build");

    assert_eq!(
        scalar
            .as_function_field_map()
            .expect("function-field pullback should build"),
        ShortWeierstrassFunctionFieldMap::new(curve.clone(), curve, field.x(), field.y(),)
            .expect("identity pullback should validate")
    );
}

#[test]
fn direct_function_field_map_matches_point_evaluation_away_from_the_kernel() {
    let curve = curve();
    let scalar = ScalarMultiplicationIsogeny::new(curve.clone(), 2)
        .expect("scalar multiplication should build");
    let map = scalar
        .as_function_field_map()
        .expect("function-field pullback should build");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !scalar.kernel_points().contains(point))
        .expect("sample curve should have a point outside the kernel");
    let image = scalar
        .evaluate(&point)
        .expect("sample point should evaluate under [2]");

    let x_value = evaluate_short_weierstrass_function_at_point(map.x_pullback(), &point)
        .expect("non-kernel point should avoid poles in x pullback");
    let y_value = evaluate_short_weierstrass_function_at_point(map.y_pullback(), &point)
        .expect("non-kernel point should avoid poles in y pullback");

    assert_eq!(Some(&x_value), AffinePoint::x_coordinate(&image));
    assert_eq!(Some(&y_value), AffinePoint::y_coordinate(&image));
}

#[test]
fn direct_p_pullback_matches_point_evaluation_away_from_the_kernel() {
    let curve = curve();
    let scalar = ScalarMultiplicationIsogeny::new(curve.clone(), 41)
        .expect("scalar multiplication should build");
    let map = scalar
        .as_function_field_map()
        .expect("direct [p]^* should build");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !scalar.kernel_points().contains(point))
        .expect("sample curve should have a point outside the kernel");
    let image = scalar
        .evaluate(&point)
        .expect("sample point should evaluate under [p]");

    let x_value = evaluate_short_weierstrass_function_at_point(map.x_pullback(), &point)
        .expect("non-kernel point should avoid poles in x pullback");
    let y_value = evaluate_short_weierstrass_function_at_point(map.y_pullback(), &point)
        .expect("non-kernel point should avoid poles in y pullback");

    assert_eq!(Some(&x_value), AffinePoint::x_coordinate(&image));
    assert_eq!(Some(&y_value), AffinePoint::y_coordinate(&image));
}
