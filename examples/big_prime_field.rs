use elliptic_algorithms_lab::fields::{BigPrimeField, traits::AmbientField};
use num_bigint::BigUint;
use num_traits::One;

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

fn print_small_field_walkthrough() -> Result<(), Box<dyn std::error::Error>> {
    let field = BigPrimeField::new(BigUint::from(17u8))?;
    let x = field.elem_from_u64(5);
    let y = field.elem_from_u64(11);
    let sum = AmbientField::add(&field, &x, &y)?;
    let product = AmbientField::mul(&field, &x, &y)?;
    let inverse_of_x = AmbientField::inverse(&field, &x)?;
    let square = field.elem_from_u64(4);
    let square_root = field.sqrt(&square)?;
    let character = field.quadratic_character_of(&square)?;

    println!("Large prime field educational walkthrough");
    println!("======================================================");
    println!();
    println!("Small sanity-check field");
    println!("------------------------");
    println!("field            = {}", field);
    println!("modulus bits     = {}", field.modulus_bits());
    println!("x                = {}", x);
    println!("y                = {}", y);
    println!("x + y            = {}", sum);
    println!("x · y            = {}", product);
    println!("x⁻¹              = {}", inverse_of_x);
    println!(
        "x · x⁻¹          = {}",
        AmbientField::mul(&field, &x, &inverse_of_x)?
    );
    println!("quadratic sample = {}", square);
    println!("χ(4)             = {:?}", character);
    println!(
        "sqrt(4)          = {}",
        square_root
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_else(|| "none".to_string())
    );
    println!();

    Ok(())
}

fn print_large_field_summary() -> Result<(), Box<dyn std::error::Error>> {
    let curve25519 = BigPrimeField::new(curve25519_prime())?;
    let secp256k1 = BigPrimeField::new(secp256k1_prime())?;

    let curve25519_sample = curve25519.elem_from_u64(4);
    let secp_sample = secp256k1.elem_from_u64(9);

    println!("Runtime large-prime fields");
    println!("--------------------------");
    println!("curve25519 field  = {}", curve25519);
    println!("curve25519 bits   = {}", curve25519.modulus_bits());
    println!(
        "χ(4) in curve25519 = {:?}",
        curve25519.quadratic_character_of(&curve25519_sample)?
    );
    println!(
        "sqrt(4) in curve25519 = {}",
        curve25519
            .sqrt(&curve25519_sample)?
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_else(|| "none".to_string())
    );
    println!();
    println!("secp256k1 field   = {}", secp256k1);
    println!("secp256k1 bits    = {}", secp256k1.modulus_bits());
    println!(
        "χ(9) in secp256k1 = {:?}",
        secp256k1.quadratic_character_of(&secp_sample)?
    );
    println!(
        "sqrt(9) in secp256k1 = {}",
        secp256k1
            .sqrt(&secp_sample)?
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_else(|| "none".to_string())
    );
    println!();
    println!(
        "note: this example keeps the story educational: explicit runtime fields, canonical residues, and exact arithmetic over named large primes."
    );

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    print_small_field_walkthrough()?;
    print_large_field_summary()?;
    Ok(())
}
