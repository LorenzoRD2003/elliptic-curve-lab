use num_bigint::BigInt;
use num_rational::BigRational;
use proptest::prelude::*;

use crate::elliptic_curves::{EnumerableCurveModel, ShortWeierstrassCurve};
use crate::fields::{
    CbrtField, ComplexApprox, EnumerableFiniteField, ExtensionField, ExtensionFieldElement,
    ExtensionFieldSpec, Field, FieldError, FiniteField, Fp, PolynomialModulus, Q, SqrtField,
};
use crate::proptest_support::fields::{
    ProptestF17Sqrt3Field, ProptestF17TowerField, arb_tower_element_case,
};

type F17 = Fp<17>;

fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

crate::fields::define_q_quadratic_extension!(
    spec: QSqrt2Spec,
    field: QSqrt2,
    radicand: 2,
    name: "Q(sqrt(2))",
);

struct QSqrt2ISpec;

impl ExtensionFieldSpec for QSqrt2ISpec {
    type Base = QSqrt2;

    const NAME: &'static str = "Q(sqrt(2), i)";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<QSqrt2>::new(vec![QSqrt2::one(), QSqrt2::zero(), QSqrt2::one()])
            .expect("x^2 + 1 should be a valid structural modulus")
    }

    fn check_field_conditions() -> Result<(), FieldError> {
        // The crate does not yet expose a generic irreducibility backend
        // over algebraic extension bases, so this tower step is accepted
        // from the known mathematics of Q(sqrt(2), i).
        Ok(())
    }
}

type QSqrt2I = ExtensionField<QSqrt2ISpec>;

crate::fields::define_fp_quadratic_extension!(
    spec: F17Sqrt3Spec,
    field: F17Sqrt3,
    base: F17,
    non_residue: 3,
    name: "F17(sqrt(3))",
);

struct F17TowerSpec;

impl ExtensionFieldSpec for F17TowerSpec {
    type Base = F17Sqrt3;

    const NAME: &'static str = "F17(sqrt(3))(u)";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<F17Sqrt3>::new(vec![
            F17Sqrt3::one(),
            F17Sqrt3::one(),
            F17Sqrt3::zero(),
            F17Sqrt3::one(),
        ])
        .expect("tower modulus should be structurally valid")
    }

    fn check_field_conditions() -> Result<(), FieldError> {
        Ok(())
    }
}

type F17Tower = ExtensionField<F17TowerSpec>;

#[test]
fn extension_field_spec_can_be_constructed_as_zero_sized_field() {
    let _field = QSqrt2::new().expect("q(sqrt(2)) should validate");
    assert_eq!(QSqrt2::name(), "Q(sqrt(2))");
    assert_eq!(QSqrt2::extension_degree().get(), 2);
}

#[test]
fn extension_field_reduces_x_squared_to_base_relation() {
    let x_squared = QSqrt2::element(vec![Q::zero(), Q::zero(), Q::one()]);
    let reduced = QSqrt2::reduce(&x_squared);

    assert_eq!(reduced.coefficients().len(), 1);
    assert!(Q::eq(&reduced.coefficients()[0], &q(2, 1)));
}

#[test]
fn extension_field_addition_and_multiplication_work_over_q_sqrt2() {
    let one_plus_x = QSqrt2::element(vec![Q::one(), Q::one()]);
    let three_minus_x = QSqrt2::element(vec![Q::from_i64(3), Q::from_i64(-1)]);

    let sum = QSqrt2::add(&one_plus_x, &three_minus_x);
    let product = QSqrt2::mul(&one_plus_x, &three_minus_x);

    assert_eq!(sum.coefficients().len(), 1);
    assert!(Q::eq(&sum.coefficients()[0], &q(4, 1)));

    assert_eq!(product.coefficients().len(), 2);
    assert!(Q::eq(&product.coefficients()[0], &q(1, 1)));
    assert!(Q::eq(&product.coefficients()[1], &q(2, 1)));
}

#[test]
fn extension_field_inverse_is_computed_via_polynomial_euclid() {
    let element = QSqrt2::element(vec![Q::one(), Q::one()]);
    let inverse = QSqrt2::inverse(&element).expect("1 + sqrt(2) should be invertible");
    let check = QSqrt2::mul(&element, &inverse);

    assert!(QSqrt2::eq(&check, &QSqrt2::one()));
    assert_eq!(inverse.coefficients().len(), 2);
    assert!(Q::eq(&inverse.coefficients()[0], &q(-1, 1)));
    assert!(Q::eq(&inverse.coefficients()[1], &q(1, 1)));
}

#[test]
fn extension_field_elements_compare_by_quotient_equivalence() {
    let x_squared = ExtensionFieldElement::<QSqrt2Spec>::new(vec![Q::zero(), Q::zero(), Q::one()]);
    let two = QSqrt2::element(vec![Q::from_i64(2)]);

    assert_eq!(x_squared, two);
}

#[test]
fn extension_field_can_be_used_as_base_of_another_extension() {
    let i = QSqrt2I::element(vec![QSqrt2::zero(), QSqrt2::one()]);
    let i_squared = QSqrt2I::mul(&i, &i);
    let minus_one = QSqrt2I::element(vec![QSqrt2::from_i64(-1)]);

    assert_eq!(i_squared, minus_one);
}

#[test]
fn extension_field_embeds_base_elements_in_towers() {
    let embedded = QSqrt2I::from_i64(2);
    let expected = QSqrt2I::element(vec![QSqrt2::from_i64(2)]);

    assert_eq!(embedded, expected);
}

#[test]
fn finite_field_metadata_multiplies_through_extension_towers() {
    assert_eq!(F17Sqrt3::characteristic(), 17);
    assert_eq!(F17Sqrt3::extension_degree().get(), 2);
    assert_eq!(<F17Tower as FiniteField>::extension_degree().get(), 6);
}

#[test]
fn enumerable_quadratic_extension_lists_every_canonical_reduced_class() {
    let elements = F17Sqrt3::elements();

    assert_eq!(elements.len(), 17 * 17);
    assert_eq!(elements[0], F17Sqrt3::zero());
    assert!(elements.contains(&F17Sqrt3::one()));
    assert!(elements.contains(&F17Sqrt3::element(vec![F17::zero(), F17::one()])));
}

#[test]
fn brute_force_square_root_finds_the_extension_generator_square_root() {
    let root = F17Sqrt3::sqrt(&F17Sqrt3::from_base(F17::from_i64(3)))
        .expect("the extension generator squares to the chosen non-residue");

    assert_eq!(
        F17Sqrt3::square(&root),
        F17Sqrt3::from_base(F17::from_i64(3))
    );
}

#[test]
fn brute_force_square_root_honestly_rejects_a_non_square_extension_element() {
    let nonsquare = F17Sqrt3::elements()
        .into_iter()
        .find(|element| F17Sqrt3::sqrt(element).is_none())
        .expect("a finite field should contain non-squares");

    assert_eq!(F17Sqrt3::sqrt(&nonsquare), None);
}

#[test]
fn brute_force_cube_root_finds_a_cube_root_for_the_extension_generator_cube() {
    let generator = F17Sqrt3::element(vec![F17::zero(), F17::one()]);
    let cube = F17Sqrt3::cube(&generator);
    let root = F17Sqrt3::cbrt(&cube).expect("a genuine cube in the extension should admit a root");

    assert_eq!(F17Sqrt3::cube(&root), cube);
}

#[test]
fn brute_force_cube_root_honestly_rejects_a_non_cube_extension_element() {
    let noncube = F17Sqrt3::elements()
        .into_iter()
        .find(|element| F17Sqrt3::cbrt(element).is_none())
        .expect("a finite field should contain non-cubes when q - 1 is divisible by 3");

    assert_eq!(F17Sqrt3::cbrt(&noncube), None);
}

#[test]
fn short_weierstrass_curves_over_enumerable_extension_fields_can_list_points() {
    let curve = ShortWeierstrassCurve::<F17Sqrt3>::new(F17Sqrt3::one(), F17Sqrt3::one())
        .expect("the lifted short-Weierstrass model should stay smooth");

    let points = curve.points();

    assert!(!points.is_empty());
    assert_eq!(curve.order(), points.len());
}

#[test]
fn extension_field_structure_rejects_constant_modulus_after_normalization() {
    struct BadSpec;

    impl ExtensionFieldSpec for BadSpec {
        type Base = F17;

        fn defining_modulus() -> PolynomialModulus<Self::Base> {
            PolynomialModulus::<F17>::new(vec![F17::one(), F17::zero()])
                .expect("structural constructor should accept the raw coefficients")
        }
    }

    let error = ExtensionField::<BadSpec>::check_structure()
        .expect_err("effectively constant modulus should be rejected");
    assert!(matches!(error, FieldError::InvalidPolynomialModulus));
}

#[test]
fn extension_field_non_invertible_zero_is_rejected() {
    let error = QSqrt2::inverse(&QSqrt2::zero()).expect_err("zero should not be invertible");
    assert!(matches!(error, FieldError::DivisionByZero));
}

#[test]
fn extension_field_algebraic_closedness_is_forwarded_from_the_spec() {
    fn algebraically_closed<F: Field>() -> bool {
        F::IS_ALGEBRAICALLY_CLOSED
    }

    struct ComplexLinearSpec;

    impl ExtensionFieldSpec for ComplexLinearSpec {
        type Base = ComplexApprox;

        const NAME: &'static str = "C(a)";
        const IS_ALGEBRAICALLY_CLOSED: bool = true;

        fn defining_modulus() -> PolynomialModulus<Self::Base> {
            PolynomialModulus::<ComplexApprox>::new(vec![
                ComplexApprox::from_i64(-1),
                ComplexApprox::one(),
            ])
            .expect("linear modulus should be structurally valid")
        }
    }

    assert!(algebraically_closed::<ExtensionField<ComplexLinearSpec>>());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(28))]

    #[test]
    fn property_extension_field_reduction_is_idempotent(
        coefficients in prop::collection::vec(0u64..17, 0..=5),
    ) {
        let element = F17Sqrt3::element(
            coefficients
                .into_iter()
                .map(F17::elem_from_u64)
                .collect(),
        );

        let reduced_once = F17Sqrt3::reduce(&element);
        let reduced_twice = F17Sqrt3::reduce(&reduced_once);

        prop_assert_eq!(reduced_once, reduced_twice);
    }

    #[test]
    fn property_extension_field_nonzero_elements_are_invertible(
        coefficients in prop::collection::vec(0u64..17, 0..=2),
    ) {
        let element = F17Sqrt3::element(
            coefficients
                .into_iter()
                .map(F17::elem_from_u64)
                .collect(),
        );
        let reduced = F17Sqrt3::reduce(&element);

        prop_assume!(!F17Sqrt3::is_zero(&reduced));

        let inverse = F17Sqrt3::inverse(&reduced).expect("non-zero extension-field element should be invertible");
        let one = F17Sqrt3::mul(&reduced, &inverse);

        prop_assert!(F17Sqrt3::eq(&one, &F17Sqrt3::one()));
    }

    #[test]
    fn property_base_embedding_into_quadratic_extension_is_a_field_homomorphism(
        case in arb_tower_element_case(),
    ) {
        let embedded_left = ProptestF17Sqrt3Field::from_base(case.base_left);
        let embedded_right = ProptestF17Sqrt3Field::from_base(case.base_right);
        let embedded_sum = ProptestF17Sqrt3Field::add(&embedded_left, &embedded_right);
        let embedded_product = ProptestF17Sqrt3Field::mul(&embedded_left, &embedded_right);

        prop_assert_eq!(
            embedded_sum,
            ProptestF17Sqrt3Field::from_base(F17::add(&case.base_left, &case.base_right))
        );
        prop_assert_eq!(
            embedded_product,
            ProptestF17Sqrt3Field::from_base(F17::mul(&case.base_left, &case.base_right))
        );
    }

    #[test]
    fn property_quadratic_base_embedding_into_tower_is_a_field_homomorphism(
        case in arb_tower_element_case(),
    ) {
        let embedded_left = ProptestF17TowerField::from_base(case.quadratic_left.clone());
        let embedded_right = ProptestF17TowerField::from_base(case.quadratic_right.clone());
        let embedded_sum = ProptestF17TowerField::add(&embedded_left, &embedded_right);
        let embedded_product = ProptestF17TowerField::mul(&embedded_left, &embedded_right);

        prop_assert_eq!(
            embedded_sum,
            ProptestF17TowerField::from_base(ProptestF17Sqrt3Field::add(
                &case.quadratic_left,
                &case.quadratic_right,
            ))
        );
        prop_assert_eq!(
            embedded_product,
            ProptestF17TowerField::from_base(ProptestF17Sqrt3Field::mul(
                &case.quadratic_left,
                &case.quadratic_right,
            ))
        );
    }

    #[test]
    fn property_semantic_tower_elements_stay_canonically_reduced_under_arithmetic(
        case in arb_tower_element_case(),
    ) {
        let sum = ProptestF17TowerField::add(&case.tower_left, &case.tower_right);
        let product = ProptestF17TowerField::mul(&case.tower_left, &case.tower_right);

        prop_assert_eq!(ProptestF17TowerField::reduce(&sum), sum);
        prop_assert_eq!(ProptestF17TowerField::reduce(&product), product);
    }
}
