use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::BigPrimeField;
use crate::fields::{
    FieldError,
    big_prime_field::BigPrimeFieldElem,
    traits::{AmbientField, QuadraticCharacterValue},
};

impl BigPrimeField {
    /// Returns one square root of `value` in `F_p` when it exists.
    ///
    /// Behavior by case:
    ///
    /// - for `p = 2`, every element is its own square root
    /// - for odd primes, this uses the Tonelli-Shanks algorithm
    ///
    /// This is the large-prime runtime analogue of the existing `Fp<P>`
    /// square-root route. The current non-residue search still scans small
    /// integers `2, 3, 4, …` until it finds a witness.
    pub fn sqrt(&self, value: &BigPrimeFieldElem) -> Result<Option<BigPrimeFieldElem>, FieldError> {
        if let Some(immediate_root) = self.sqrt_immediate_case(value) {
            return Ok(Some(immediate_root));
        }

        if self.quadratic_character_of(value)? == QuadraticCharacterValue::NonResidue {
            return Ok(None);
        }

        if self.modulus_congruent_to_three_mod_four() {
            return Ok(Some(self.sqrt_when_modulus_is_three_mod_four(value)));
        }

        self.sqrt_via_tonelli_shanks(value).map(Some)
    }

    /// Returns one square root together with its additive inverse.
    pub fn sqrt_pair(
        &self,
        value: &BigPrimeFieldElem,
    ) -> Result<Option<(BigPrimeFieldElem, BigPrimeFieldElem)>, FieldError> {
        let root = self.sqrt(value)?;
        Ok(root.map(|left| {
            let right = AmbientField::neg(self, &left);
            (left, right)
        }))
    }

    fn sqrt_immediate_case(&self, value: &BigPrimeFieldElem) -> Option<BigPrimeFieldElem> {
        if AmbientField::is_zero(self, value) {
            return Some(AmbientField::zero(self));
        }

        if self.modulus() == &BigUint::from(2u8) {
            return Some(value.clone());
        }
        None
    }

    fn modulus_congruent_to_three_mod_four(&self) -> bool {
        self.modulus() % BigUint::from(4u8) == BigUint::from(3u8)
    }

    fn sqrt_when_modulus_is_three_mod_four(&self, value: &BigPrimeFieldElem) -> BigPrimeFieldElem {
        let exponent = (self.modulus() + BigUint::one()) >> 2usize;
        self.pow_biguint(value, &exponent)
    }

    fn sqrt_via_tonelli_shanks(
        &self,
        value: &BigPrimeFieldElem,
    ) -> Result<BigPrimeFieldElem, FieldError> {
        let (odd_part, two_adic_exponent) = decompose_modulus_minus_one(self.modulus());
        let non_residue = self.find_quadratic_non_residue()?;
        let mut m = two_adic_exponent;
        let mut c = self.pow_biguint(&non_residue, &odd_part);
        let mut t = self.pow_biguint(value, &odd_part);
        let half_odd_part_plus_one = (&odd_part + BigUint::one()) >> 1usize;
        let mut r = self.pow_biguint(value, &half_odd_part_plus_one);

        while !AmbientField::eq(self, &t, &AmbientField::one(self)) {
            let Some(next_index) = self.first_power_of_t_landing_at_one(&t, m) else {
                return Err(FieldError::Unsupported(
                    "Tonelli-Shanks failed to find the expected 2-power descent index in the current big prime-field backend",
                ));
            };

            let tonelli_step = self.tonelli_shanks_step(&c, &r, &t, m, next_index);
            m = next_index;
            c = tonelli_step.next_c;
            r = tonelli_step.next_r;
            t = tonelli_step.next_t;
        }

        Ok(r)
    }

    fn find_quadratic_non_residue(&self) -> Result<BigPrimeFieldElem, FieldError> {
        let mut candidate = BigUint::from(2u8);

        loop {
            let candidate_mod_p = self.elem(candidate.clone());
            if self.quadratic_character_of(&candidate_mod_p)? == QuadraticCharacterValue::NonResidue
            {
                return Ok(candidate_mod_p);
            }

            candidate += BigUint::one();
        }
    }

    fn first_power_of_t_landing_at_one(
        &self,
        t: &BigPrimeFieldElem,
        exponent_bound: u32,
    ) -> Option<u32> {
        let mut i = 1u32;
        let mut t_power =
            AmbientField::mul(self, t, t).expect("prime-field multiplication should stay total");

        while i < exponent_bound && !AmbientField::eq(self, &t_power, &AmbientField::one(self)) {
            t_power = AmbientField::mul(self, &t_power, &t_power)
                .expect("prime-field multiplication should stay total");
            i += 1;
        }

        (i < exponent_bound).then_some(i)
    }

    fn tonelli_shanks_step(
        &self,
        c: &BigPrimeFieldElem,
        r: &BigPrimeFieldElem,
        t: &BigPrimeFieldElem,
        current_exponent: u32,
        next_index: u32,
    ) -> TonelliShanksStep {
        let exponent = BigUint::one() << (current_exponent - next_index - 1);
        let b = self.pow_biguint(c, &exponent);
        let next_c =
            AmbientField::mul(self, &b, &b).expect("prime-field multiplication should stay total");
        let next_r =
            AmbientField::mul(self, r, &b).expect("prime-field multiplication should stay total");
        let next_t = AmbientField::mul(self, t, &next_c)
            .expect("prime-field multiplication should stay total");

        TonelliShanksStep {
            next_c,
            next_r,
            next_t,
        }
    }
}

struct TonelliShanksStep {
    next_c: BigPrimeFieldElem,
    next_r: BigPrimeFieldElem,
    next_t: BigPrimeFieldElem,
}

fn decompose_modulus_minus_one(modulus: &BigUint) -> (BigUint, u32) {
    let mut odd_part = modulus - BigUint::one();
    let mut two_adic_exponent = 0u32;

    while (&odd_part & BigUint::one()).is_zero() {
        odd_part >>= 1usize;
        two_adic_exponent += 1;
    }

    (odd_part, two_adic_exponent)
}
