use num_bigint::BigUint;
use num_traits::One;

use crate::fields::{
    FieldError, Fp,
    big_prime_field::BigPrimeField,
    traits::{
        AmbientField, Field, QuadraticCharacterFiniteField, QuadraticCharacterValue, SqrtField,
    },
};

type F17 = Fp<17>;
type F2 = Fp<2>;

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn curve25519_prime() -> BigUint {
    BigUint::from(2u8).pow(255) - BigUint::from(19u8)
}

fn secp256k1_prime() -> BigUint {
    BigUint::from(2u8).pow(256)
        - BigUint::from(2u8).pow(32)
        - BigUint::from(2u16).pow(9)
        - BigUint::from(2u8).pow(8)
        - BigUint::from(2u8).pow(7)
        - BigUint::from(2u8).pow(6)
        - BigUint::from(2u8).pow(4)
        - BigUint::one()
}

#[test]
fn rejects_moduli_below_two() {
    assert_eq!(
        BigPrimeField::new(BigUint::from(0u8)),
        Err(FieldError::InvalidBigModulus {
            modulus: "0".into()
        })
    );
    assert_eq!(
        BigPrimeField::new(BigUint::from(1u8)),
        Err(FieldError::InvalidBigModulus {
            modulus: "1".into()
        })
    );
}

#[test]
fn rejects_composite_runtime_modulus() {
    assert_eq!(
        BigPrimeField::new(bu(21)),
        Err(FieldError::InvalidBigModulus {
            modulus: "21".into()
        })
    );
}

#[test]
fn accepts_small_and_large_prime_moduli() {
    let small = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let large =
        BigPrimeField::new(curve25519_prime()).expect("2^255 - 19 should define a prime field");
    let secp = BigPrimeField::new(secp256k1_prime()).expect("the secp256k1 prime should be valid");

    assert_eq!(small.modulus(), &bu(17));
    assert!(large.modulus().bits() >= 255);
    assert!(secp.modulus().bits() >= 256);
}

#[test]
fn display_and_bit_length_expose_a_small_educational_summary() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let curve25519 =
        BigPrimeField::new(curve25519_prime()).expect("curve25519 prime should be valid");

    assert_eq!(field.to_string(), "F_17");
    assert_eq!(field.modulus_bits(), 5);
    assert!(curve25519.to_string().starts_with("F_"));
    assert!(curve25519.modulus_bits() >= 255);
}

#[test]
fn canonicalizes_elements_by_reduction_mod_p() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let reduced = field.elem(bu(52));

    assert_eq!(reduced.value(), &bu(1));
    assert!(field.is_canonical_value(reduced.value()));
    assert!(!field.is_canonical_value(&bu(17)));
}

#[test]
fn zero_and_one_use_the_canonical_representatives() {
    let field = BigPrimeField::new(curve25519_prime()).expect("prime field should be valid");

    assert_eq!(AmbientField::zero(&field).value(), &BigUint::from(0u8));
    assert_eq!(AmbientField::one(&field).value(), &BigUint::from(1u8));
}

#[test]
fn equal_residue_classes_share_the_same_canonical_value() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let left = field.elem(bu(35));
    let right = field.elem(bu(1));

    assert_eq!(left, right);
    assert_eq!(left.value(), &bu(1));
    assert!(AmbientField::eq(&field, &left, &right));
}

#[test]
fn addition_matches_the_small_static_prime_field_route() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let left = field.elem_from_u64(13);
    let right = field.elem_from_u64(9);
    let sum = AmbientField::add(&field, &left, &right).expect("field addition should succeed");

    let static_sum = F17::add(&F17::elem_from_u64(13), &F17::elem_from_u64(9));

    assert_eq!(sum.value(), &BigUint::from(static_sum.value()));
}

#[test]
fn subtraction_and_negation_wrap_canonically() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let left = field.elem_from_u64(3);
    let right = field.elem_from_u64(5);

    let neg = AmbientField::neg(&field, &right);
    let difference =
        AmbientField::sub(&field, &left, &right).expect("field subtraction should succeed");

    assert_eq!(neg.value(), &bu(12));
    assert_eq!(difference.value(), &bu(15));
}

#[test]
fn multiplication_matches_the_small_static_prime_field_route() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let left = field.elem_from_u64(7);
    let right = field.elem_from_u64(11);
    let product =
        AmbientField::mul(&field, &left, &right).expect("field multiplication should succeed");

    let static_product = F17::mul(&F17::elem_from_u64(7), &F17::elem_from_u64(11));

    assert_eq!(product.value(), &BigUint::from(static_product.value()));
}

#[test]
fn inverse_and_division_recover_the_expected_unit_relations() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let numerator = field.elem_from_u64(7);
    let denominator = field.elem_from_u64(5);

    let inverse =
        AmbientField::inverse(&field, &denominator).expect("non-zero residue should be invertible");
    let quotient = AmbientField::div(&field, &numerator, &denominator)
        .expect("division by a unit should succeed");

    assert_eq!(inverse.value(), &bu(7));
    assert_eq!(
        AmbientField::mul(&field, &denominator, &inverse).unwrap(),
        AmbientField::one(&field)
    );
    assert_eq!(
        AmbientField::mul(&field, &quotient, &denominator).unwrap(),
        numerator
    );
}

#[test]
fn inverse_and_division_reject_zero_honestly() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let one = AmbientField::one(&field);
    let zero = AmbientField::zero(&field);

    assert_eq!(
        AmbientField::inverse(&field, &zero),
        Err(FieldError::DivisionByZero)
    );
    assert_eq!(
        AmbientField::div(&field, &one, &zero),
        Err(FieldError::DivisionByZero)
    );
}

#[test]
fn ambient_field_impl_matches_the_inherent_runtime_arithmetic() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let left = field.elem_from_u64(4);
    let right = field.elem_from_u64(15);

    let trait_sum = AmbientField::add(&field, &left, &right).expect("ambient addition should work");
    let default_difference = AmbientField::sub(&field, &trait_sum, &right)
        .expect("ambient subtraction should reuse addition plus negation");

    assert_eq!(default_difference, left);
}

#[test]
fn biguint_power_matches_the_small_static_prime_field_route() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let value = field.elem_from_u64(5);
    let exponent = BigUint::from(13u8);

    let powered = field.pow_biguint(&value, &exponent);
    let static_power = F17::pow(&F17::elem_from_u64(5), 13);

    assert_eq!(powered.value(), &BigUint::from(static_power.value()));
}

#[test]
fn quadratic_character_matches_the_small_static_prime_field_route() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");

    for raw_value in 0..17 {
        let runtime_value = field.elem_from_u64(raw_value);
        let static_value = F17::elem_from_u64(raw_value);

        assert_eq!(
            field.quadratic_character_of(&runtime_value),
            F17::quadratic_character_of(&static_value)
        );
    }
}

#[test]
fn quadratic_character_reports_characteristic_two_as_unsupported() {
    let field = BigPrimeField::new(bu(2)).expect("2 should define a prime field");
    let one = field.elem_from_u64(1);

    assert_eq!(
        field.quadratic_character_of(&one),
        Err(FieldError::Unsupported(
            "quadratic character is only implemented for finite fields of odd characteristic"
        ))
    );
}

#[test]
fn has_square_root_matches_quadratic_character_on_f17() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");

    for raw_value in 0..17 {
        let value = field.elem_from_u64(raw_value);
        let character = field
            .quadratic_character_of(&value)
            .expect("F_17 should support the quadratic character");
        let has_root = field
            .has_square_root(&value)
            .expect("F_17 should support square-root existence checks");

        assert_eq!(
            has_root,
            matches!(
                character,
                QuadraticCharacterValue::Zero | QuadraticCharacterValue::Residue
            )
        );
    }
}

#[test]
fn sqrt_matches_the_small_static_prime_field_route_on_f17() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");

    for raw_value in 0..17 {
        let runtime_value = field.elem_from_u64(raw_value);
        let static_value = F17::elem_from_u64(raw_value);

        let runtime_root = field
            .sqrt(&runtime_value)
            .expect("F_17 sqrt should be supported");
        let static_root = F17::sqrt(&static_value).map(|root| BigUint::from(root.value()));

        assert_eq!(
            runtime_root.as_ref().map(|root| root.value().clone()),
            static_root
        );
    }
}

#[test]
fn sqrt_pair_returns_opposite_runtime_roots() {
    let field = BigPrimeField::new(bu(17)).expect("17 should define a prime field");
    let four = field.elem_from_u64(4);
    let (left, right) = field
        .sqrt_pair(&four)
        .expect("sqrt_pair should be supported over F_17")
        .expect("4 should be a square over F_17");

    assert_eq!(AmbientField::mul(&field, &left, &left).unwrap(), four);
    assert_eq!(AmbientField::mul(&field, &right, &right).unwrap(), four);
    assert_eq!(right, AmbientField::neg(&field, &left));
}

#[test]
fn sqrt_in_characteristic_two_returns_the_input() {
    let field = BigPrimeField::new(bu(2)).expect("2 should define a prime field");
    let one = field.elem_from_u64(1);

    assert_eq!(field.sqrt(&one), Ok(Some(one.clone())));
    assert_eq!(field.sqrt_pair(&one), Ok(Some((one.clone(), one))));
    assert_eq!(F2::sqrt(&F2::one()).map(|root| root.value()), Some(1));
}
