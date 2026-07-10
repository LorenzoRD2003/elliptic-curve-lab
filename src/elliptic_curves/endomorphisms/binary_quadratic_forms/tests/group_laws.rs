use crate::elliptic_curves::endomorphisms::{
    BinaryQuadraticForm, QuadraticClassGroup, quadratic_orders::QuadraticDiscriminant,
};

use super::z;

#[test]
fn public_compose_satisfies_group_laws_for_small_discriminants() {
    for discriminant in [-3, -4, -20, -23, -31, -84] {
        assert_group_laws_for(discriminant);
    }
}

fn assert_group_laws_for(discriminant: i64) {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(discriminant))
        .expect("test discriminants should define imaginary quadratic orders");
    let forms = class_group.enumerate_reduced_forms();
    let identity = principal_form(discriminant);

    assert!(forms.contains(&identity));

    for left in &forms {
        for right in &forms {
            let product = class_group
                .compose(left, right)
                .expect("enumerated class-group representatives should compose");

            assert!(
                forms.contains(&product),
                "product {product:?} for D = {discriminant} should be enumerated"
            );
        }
    }

    for form in &forms {
        assert_eq!(
            class_group
                .compose(&identity, form)
                .expect("identity should compose on the left"),
            *form
        );
        assert_eq!(
            class_group
                .compose(form, &identity)
                .expect("identity should compose on the right"),
            *form
        );

        let inverse = form
            .conjugate()
            .reduce_positive_definite()
            .expect("positive-definite conjugates should reduce");

        assert!(forms.contains(&inverse));
        assert_eq!(
            class_group
                .compose(form, &inverse)
                .expect("form should compose with its conjugate"),
            identity
        );
        assert_eq!(
            class_group
                .compose(&inverse, form)
                .expect("conjugate should compose with the form"),
            identity
        );
    }

    for left in &forms {
        for middle in &forms {
            for right in &forms {
                let left_associated = class_group
                    .compose(
                        &class_group
                            .compose(left, middle)
                            .expect("left product should compose"),
                        right,
                    )
                    .expect("left-associated product should compose");
                let right_associated = class_group
                    .compose(
                        left,
                        &class_group
                            .compose(middle, right)
                            .expect("right product should compose"),
                    )
                    .expect("right-associated product should compose");

                assert_eq!(left_associated, right_associated);
            }
        }
    }
}

fn principal_form(discriminant: i64) -> BinaryQuadraticForm {
    if discriminant.rem_euclid(4) == 0 {
        BinaryQuadraticForm::new(z(1), z(0), z(-discriminant / 4))
    } else {
        BinaryQuadraticForm::new(z(1), z(1), z((1 - discriminant) / 4))
    }
}
