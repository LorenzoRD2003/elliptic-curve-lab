//! Spike tests for future static large prime fields.
//!
//! These tests intentionally exercise `crypto-bigint` directly instead of
//! defining a new public field backend. The milestone question is whether
//! constant Montgomery arithmetic over `ℤ/pℤ` gives us the right raw material
//! before we widen `Field`/`FiniteField` metadata.

use crypto_bigint::{
    U256, const_prime_monty_params,
    modular::{ConstMontyForm, ConstPrimeMontyParams},
};

const P256_MODULUS_HEX: &str = "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff";

const_prime_monty_params!(
    P256Prime,
    U256,
    P256_MODULUS_HEX,
    6,
    "P-256 prime modulus used only by the repo-local crypto-bigint spike."
);

type P256Elem = ConstMontyForm<P256Prime, { U256::LIMBS }>;

fn elem(value: u64) -> P256Elem {
    P256Elem::new(&U256::from_u64(value))
}

fn repr(element: &P256Elem) -> U256 {
    element.retrieve()
}

#[test]
fn const_monty_form_supports_basic_prime_field_arithmetic() {
    let three = elem(3);
    let five = elem(5);

    assert_eq!(repr(&(&three + &five)), U256::from_u64(8));
    assert_eq!(repr(&(&three * &five)), U256::from_u64(15));
    assert_eq!(
        repr(&(-&three)),
        U256::from_be_hex(P256_MODULUS_HEX) - U256::from_u64(3)
    );

    let inverse = five.invert().expect("5 is non-zero modulo the P-256 prime");
    assert_eq!(repr(&(&five * &inverse)), U256::ONE);

    let power = three.pow(&U256::from_u64(20));
    assert_eq!(repr(&power), U256::from_u64(3_486_784_401));
}

#[test]
fn const_prime_monty_params_unlock_square_roots() {
    let four = elem(4);
    assert_eq!(four.sqrt().expect("4 is a square modulo p"), elem(2));

    let generator = elem(u64::from(P256Prime::PRIME_PARAMS.generator().get()));
    assert!(generator.sqrt().is_none().to_bool_vartime());
}

#[test]
fn const_monty_form_retrieves_canonical_representatives() {
    let modulus = U256::from_be_hex(P256_MODULUS_HEX);

    assert_eq!(repr(&P256Elem::new(&modulus)), U256::ZERO);
    assert_eq!(
        repr(&P256Elem::new(&U256::from_u64(42))),
        U256::from_u64(42)
    );
}
