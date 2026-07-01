use num_bigint::BigUint;

use crate::numerics::quadratic_forms::{
    DiagonalBinaryQuadraticForm, DiagonalBinaryQuadraticRepresentation, QuadraticFormError,
};

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

#[test]
fn diagonal_form_rejects_zero_coefficient() {
    assert_eq!(
        DiagonalBinaryQuadraticForm::new(bu(0)),
        Err(QuadraticFormError::ZeroDiagonalCoefficient)
    );
}

#[test]
fn primitive_representations_wrap_cornacchia_solutions() {
    let form = DiagonalBinaryQuadraticForm::new(bu(1)).expect("x² + y² is positive diagonal");
    let representations = form
        .primitive_representations(&bu(65))
        .expect("65 should have primitive representations as a sum of two squares");

    assert_eq!(form.coefficient(), &bu(1));
    assert_eq!(
        representations,
        vec![
            DiagonalBinaryQuadraticRepresentation::new(bu(7), bu(4)),
            DiagonalBinaryQuadraticRepresentation::new(bu(8), bu(1)),
        ]
    );
    for representation in representations {
        assert_eq!(representation.value(form.coefficient()), bu(65));
    }
}

#[test]
fn primitive_representations_filter_non_primitive_cornacchia_candidates() {
    let form = DiagonalBinaryQuadraticForm::new(bu(4)).expect("d = 4 is positive");

    assert_eq!(
        form.primitive_representations(&bu(20)),
        Ok(vec![DiagonalBinaryQuadraticRepresentation::new(
            bu(4),
            bu(1)
        )])
    );
}

#[test]
fn primitively_represents_reports_whether_a_primitive_representation_exists() {
    let form = DiagonalBinaryQuadraticForm::new(bu(3)).expect("d = 3 is positive");

    assert!(
        form.primitively_represents(&bu(13))
            .expect("13 = 1² + 3·2²")
    );
    assert!(
        !form
            .primitively_represents(&bu(11))
            .expect("11 has no primitive representation by x² + 3y²")
    );
}
