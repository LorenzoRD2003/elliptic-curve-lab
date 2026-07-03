use crate::elliptic_curves::short_weierstrass::rational_torsion::{
    RationalTorsionError, RationalTorsionGroup, RationalTorsionGroupShape,
    enumeration::{MAZUR_CYCLIC_ORDERS, MAZUR_PRODUCT_PARAMETERS, MAZUR_TORSION_EXPONENT_BOUND},
};

#[test]
fn mazur_constants_record_the_first_route_bounds() {
    assert_eq!(MAZUR_TORSION_EXPONENT_BOUND, 27_720);
    assert_eq!(MAZUR_CYCLIC_ORDERS, &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]);
    assert_eq!(MAZUR_PRODUCT_PARAMETERS, &[1, 2, 3, 4]);
}

#[test]
fn rational_torsion_group_cardinality_matches_mazur_shape() {
    assert_eq!(
        RationalTorsionGroup::new(RationalTorsionGroupShape::Trivial)
            .expect("trivial group is a Mazur shape")
            .cardinality(),
        1
    );
    assert_eq!(
        RationalTorsionGroup::new(RationalTorsionGroupShape::Cyclic { order: 6 })
            .expect("order 6 is a Mazur cyclic shape")
            .cardinality(),
        6
    );
    assert_eq!(
        RationalTorsionGroup::new(RationalTorsionGroupShape::ProductZ2Z2m { m: 2 })
            .expect("m = 2 is a Mazur product shape")
            .cardinality(),
        8
    );
}

#[test]
fn rational_torsion_group_constructor_rejects_non_mazur_shapes() {
    assert_eq!(
        RationalTorsionGroup::new(RationalTorsionGroupShape::Cyclic { order: 11 }),
        Err(RationalTorsionError::InvalidMazurShape {
            family: "cyclic",
            value: 11,
        })
    );
    assert_eq!(
        RationalTorsionGroup::new(RationalTorsionGroupShape::ProductZ2Z2m { m: 5 }),
        Err(RationalTorsionError::InvalidMazurShape {
            family: "product",
            value: 5,
        })
    );
}
