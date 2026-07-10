use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, QuadraticClassGroup},
    quadratic_orders::QuadraticDiscriminant,
};

use super::z;

#[test]
fn gp_dirichlet_composition_tables_are_ready_for_compose_tests() {
    // External source: GP/PARI `qfbcomp(Qfb(a,b,c), Qfb(a',b',c'))`.
    //
    // These fixtures exercise the public `compose()` wrapper against small
    // GP/PARI tables.
    assert_gp_composition_table(-20, &[(1, 0, 5), (2, 2, 3)], &[&[0, 1], &[1, 0]]);
    assert_gp_composition_table(
        -23,
        &[(1, 1, 6), (2, -1, 3), (2, 1, 3)],
        &[&[0, 1, 2], &[1, 2, 0], &[2, 0, 1]],
    );
    assert_gp_composition_table(
        -31,
        &[(1, 1, 8), (2, -1, 4), (2, 1, 4)],
        &[&[0, 1, 2], &[1, 2, 0], &[2, 0, 1]],
    );
    assert_gp_composition_table(
        -84,
        &[(1, 0, 21), (2, 2, 11), (3, 0, 7), (5, 4, 5)],
        &[&[0, 1, 2, 3], &[1, 0, 3, 2], &[2, 3, 0, 1], &[3, 2, 1, 0]],
    );
}

fn assert_gp_composition_table(
    discriminant: i64,
    coefficient_rows: &[(i64, i64, i64)],
    products: &[&[usize]],
) {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(discriminant))
        .expect("test discriminants should define imaginary quadratic orders");
    let forms = coefficient_rows
        .iter()
        .map(|(a, b, c)| BinaryQuadraticForm::new(z(*a), z(*b), z(*c)))
        .collect::<Vec<_>>();

    assert_eq!(class_group.enumerate_reduced_forms(), forms);
    assert_eq!(products.len(), forms.len());

    for row in products {
        assert_eq!(row.len(), forms.len());
        assert!(row.iter().all(|&index| index < forms.len()));
    }

    for index in 0..forms.len() {
        assert_eq!(products[0][index], index);
        assert_eq!(products[index][0], index);
    }

    for (index, form) in forms.iter().enumerate() {
        let inverse = form
            .conjugate()
            .reduce_positive_definite()
            .expect("conjugates of positive-definite forms should reduce");
        let inverse_index = forms
            .iter()
            .position(|candidate| candidate == &inverse)
            .expect("the GP table should contain the reduced inverse");

        assert_eq!(products[index][inverse_index], 0);
        assert_eq!(products[inverse_index][index], 0);
    }

    for left in 0..forms.len() {
        for middle in 0..forms.len() {
            for right in 0..forms.len() {
                let left_associated = products[products[left][middle]][right];
                let right_associated = products[left][products[middle][right]];

                assert_eq!(left_associated, right_associated);
            }
        }
    }

    for (left_index, left) in forms.iter().enumerate() {
        for (right_index, right) in forms.iter().enumerate() {
            let product = class_group
                .compose(left, right)
                .expect("GP fixture forms should compose by Dirichlet/Gauss");

            assert_eq!(product, forms[products[left_index][right_index]]);
        }
    }
}
