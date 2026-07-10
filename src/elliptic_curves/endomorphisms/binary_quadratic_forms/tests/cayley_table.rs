use crate::elliptic_curves::endomorphisms::{
    BinaryQuadraticForm, QuadraticClassGroup, quadratic_orders::QuadraticDiscriminant,
};

use super::z;

#[test]
fn cayley_table_records_representatives_and_product_indices() {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-84))
        .expect("D = -84 should be supported");

    let table = class_group
        .cayley_table()
        .expect("small enumerated class group should produce a table");

    assert_eq!(table.discriminant().value(), &z(-84));
    assert_eq!(
        table.representatives(),
        &[
            BinaryQuadraticForm::new(z(1), z(0), z(21)),
            BinaryQuadraticForm::new(z(2), z(2), z(11)),
            BinaryQuadraticForm::new(z(3), z(0), z(7)),
            BinaryQuadraticForm::new(z(5), z(4), z(5)),
        ]
    );
    assert_eq!(
        table.products(),
        &[
            vec![0, 1, 2, 3],
            vec![1, 0, 3, 2],
            vec![2, 3, 0, 1],
            vec![3, 2, 1, 0],
        ]
    );
    assert_eq!(table.class_number(), 4);
    assert_eq!(
        table.representative(2),
        Some(&BinaryQuadraticForm::new(z(3), z(0), z(7)))
    );
    assert_eq!(table.product_index(1, 2), Some(3));
    assert_eq!(table.product_index(4, 0), None);
}
