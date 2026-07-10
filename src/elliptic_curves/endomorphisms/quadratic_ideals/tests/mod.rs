use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::endomorphisms::quadratic_orders::{
    ImaginaryQuadraticOrder, QuadraticDiscriminant,
};

mod errors;
mod gp_ideal_form;
mod prime_behavior;
mod prime_norm_ideal;

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn z(value: i64) -> BigInt {
    BigInt::from(value)
}

fn maximal_order(discriminant: i64) -> ImaginaryQuadraticOrder {
    let quadratic_discriminant = QuadraticDiscriminant::new(discriminant);
    ImaginaryQuadraticOrder::new(quadratic_discriminant, bu(1))
        .expect("test discriminant should define an imaginary quadratic maximal order")
}
