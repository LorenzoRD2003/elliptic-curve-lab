use proptest::prelude::*;

use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::division_polynomials::DivisionPolynomialForm,
};
use crate::fields::{Fp, traits::Field};
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_division_polynomial_case,
};

type F17 = Fp<17>;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn generated_division_polynomial_cases_match_recomputation(
        case in arb_division_polynomial_case::<17>(CurveStrategyConfig::default()),
    ) {
        let recomputed = case.curve.division_polynomial(case.index).unwrap();
        prop_assert_eq!(recomputed, case.polynomial);
    }

}

#[test]
fn odd_division_polynomials_depend_only_on_x() {
    let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
        .expect("curve should be non-singular");

    for index in [1usize, 3, 5, 7] {
        let polynomial = curve.odd_division_polynomial(index).unwrap();
        assert!(matches!(
            curve.division_polynomial(index).unwrap(),
            DivisionPolynomialForm::InX(_)
        ));
        assert_eq!(
            curve.division_polynomial(index).unwrap().x_factor(),
            &polynomial
        );
    }
}

#[test]
fn torsion_points_and_exact_order_points_are_distinguished() {
    let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
        .expect("curve should be non-singular");

    let torsion = curve.torsion_points_from_division_polynomial(6).unwrap();
    let exact = curve
        .exact_n_torsion_points_from_division_polynomial(6)
        .unwrap();

    assert!(exact.len() <= torsion.len());
    assert!(
        exact
            .iter()
            .all(|point| torsion.iter().any(|other| other == point))
    );
}
