use crate::fields::traits::*;
use num_complex::Complex64;

use crate::fields::{
    extension_field::{ExtensionField, ExtensionFieldSpec},
    polynomial_field::PolynomialModulus,
};

pub(super) type F2 = crate::fields::Fp2;
pub(super) type F3 = crate::fields::Fp3;
pub(super) type F5 = crate::fields::Fp5;
pub(super) type F7 = crate::fields::Fp7;

#[derive(Clone, Copy)]
pub(super) struct F4GeneralWeierstrassSpec;

impl ExtensionFieldSpec for F4GeneralWeierstrassSpec {
    type Base = F2;

    const NAME: &'static str = "F4 for general Weierstrass tests";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<Self::Base>::new(vec![F2::one(), F2::one(), F2::one()])
            .expect("x^2 + x + 1 should be a valid structural modulus")
    }

    fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
        Self::defining_modulus().check_field_modulus_requirements()
    }
}

pub(super) type F4 = ExtensionField<F4GeneralWeierstrassSpec>;

pub(super) fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}
