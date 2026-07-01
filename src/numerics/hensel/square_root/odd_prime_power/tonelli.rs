use num_bigint::BigUint;
use num_traits::{One, Zero};

pub(super) fn tonelli_shanks_mod_odd_prime(value_mod_p: &BigUint, p: &BigUint) -> Option<BigUint> {
    if value_mod_p.is_zero() {
        return Some(BigUint::zero());
    }

    let one = BigUint::one();
    let two = BigUint::from(2u8);
    if !is_quadratic_residue_mod_odd_prime(value_mod_p, p) {
        return None;
    }
    if p % &BigUint::from(4u8) == BigUint::from(3u8) {
        return Some(value_mod_p.modpow(&((p + &one) / &BigUint::from(4u8)), p));
    }

    let mut q = p - &one;
    let mut s = 0u32;
    while (&q % &two).is_zero() {
        q /= &two;
        s += 1;
    }

    let z = quadratic_non_residue_mod_odd_prime(p);
    let mut m = s;
    let mut c = z.modpow(&q, p);
    let mut t = value_mod_p.modpow(&q, p);
    let mut r = value_mod_p.modpow(&((&q + &one) / &two), p);

    while t != one {
        let mut i = 1u32;
        let mut t_power = (&t * &t) % p;
        while i < m && t_power != one {
            t_power = (&t_power * &t_power) % p;
            i += 1;
        }
        if i == m {
            return None;
        }

        let exponent = BigUint::one() << (m - i - 1);
        let b = c.modpow(&exponent, p);
        let b_squared = (&b * &b) % p;
        r = (&r * &b) % p;
        t = (&t * &b_squared) % p;
        c = b_squared;
        m = i;
    }

    Some(r)
}

fn is_quadratic_residue_mod_odd_prime(value_mod_p: &BigUint, p: &BigUint) -> bool {
    value_mod_p.modpow(&((p - BigUint::one()) / BigUint::from(2u8)), p) == BigUint::one()
}

fn quadratic_non_residue_mod_odd_prime(p: &BigUint) -> BigUint {
    let mut candidate = BigUint::from(2u8);
    while is_quadratic_residue_mod_odd_prime(&candidate, p) {
        candidate += BigUint::one();
    }
    candidate
}
