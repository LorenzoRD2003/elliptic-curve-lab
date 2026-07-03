use crate::elliptic_curves::short_weierstrass::rational_torsion::{
    RationalTorsionError, RationalTorsionGroup, RationalTorsionGroupShape,
    mazur::{MAZUR_CYCLIC_ORDERS, MAZUR_PRODUCT_PARAMETERS, MAZUR_TORSION_EXPONENT_BOUND},
};

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
        RationalTorsionGroup::new(RationalTorsionGroupShape::Cyclic { order: 1 }),
        Err(RationalTorsionError::InvalidMazurShape {
            family: "cyclic",
            value: 1,
        })
    );
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

#[test]
fn mazur_cyclic_orders_are_exactly_the_nontrivial_cyclic_shapes() {
    assert_eq!(MAZUR_CYCLIC_ORDERS, &[2, 3, 4, 5, 6, 7, 8, 9, 10, 12]);
    assert_eq!(MAZUR_PRODUCT_PARAMETERS, &[1, 2, 3, 4]);
    assert_eq!(MAZUR_TORSION_EXPONENT_BOUND, 27_720);
}

#[test]
fn rational_torsion_group_classifies_verified_point_orders() {
    assert_eq!(
        RationalTorsionGroup::from_verified_point_orders(&[1])
            .expect("identity-only torsion should classify")
            .shape(),
        RationalTorsionGroupShape::Trivial
    );
    assert_eq!(
        RationalTorsionGroup::from_verified_point_orders(&[1, 2, 3, 6, 6, 3])
            .expect("cyclic order-six torsion should classify")
            .shape(),
        RationalTorsionGroupShape::Cyclic { order: 6 }
    );
    assert_eq!(
        RationalTorsionGroup::from_verified_point_orders(&[1, 2, 2, 2])
            .expect("full two-torsion should classify")
            .shape(),
        RationalTorsionGroupShape::ProductZ2Z2m { m: 1 }
    );
}

#[test]
fn rational_torsion_group_rejects_impossible_verified_order_profile() {
    assert_eq!(
        RationalTorsionGroup::from_verified_point_orders(&[1, 2, 2]),
        Err(RationalTorsionError::InconsistentMazurShape { point_count: 3 })
    );
}
