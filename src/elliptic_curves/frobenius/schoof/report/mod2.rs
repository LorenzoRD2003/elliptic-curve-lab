use num_bigint::BigUint;

use crate::fields::traits::Field;
use crate::numerics::chinese_remainder::Congruence;
use crate::polynomials::DensePolynomial;

/// Report for Schoof's `ℓ = 2` trace computation.
///
/// For a short-Weierstrass curve `E: y^2 = f(x)` over an odd finite field
/// `F_q`, We can compute `t mod 2` by deciding whether `E(F_q)` contains a
/// rational point of order `2`. This is equivalent to asking whether the cubic
/// `f(x)` has a root in `F_q`, equivalently whether
///
/// `deg gcd(f(x), x^q - x) > 0`.
#[derive(Clone, Debug, PartialEq)]
pub struct SchoofTraceMod2Report<F: Field> {
    field_order: u128,
    cubic: DensePolynomial<F>,
    x_q_mod_cubic: DensePolynomial<F>,
    gcd: DensePolynomial<F>,
    has_rational_two_torsion: bool,
    trace_mod_2: u8,
}

impl<F: Field> SchoofTraceMod2Report<F> {
    pub(crate) fn new(
        field_order: u128,
        cubic: DensePolynomial<F>,
        x_q_mod_cubic: DensePolynomial<F>,
        gcd: DensePolynomial<F>,
    ) -> Self {
        let has_rational_two_torsion = gcd.degree().is_some_and(|degree| degree > 0);
        let trace_mod_2 = if has_rational_two_torsion { 0 } else { 1 };

        Self {
            field_order,
            cubic,
            x_q_mod_cubic,
            gcd,
            has_rational_two_torsion,
            trace_mod_2,
        }
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> u128 {
        self.field_order
    }

    /// Returns the short-Weierstrass cubic `f(x) = x^3 + ax + b`.
    pub fn cubic(&self) -> &DensePolynomial<F> {
        &self.cubic
    }

    /// Returns the reduced class of `x^q` in `F_q[x] / (f(x))`.
    pub fn x_q_mod_cubic(&self) -> &DensePolynomial<F> {
        &self.x_q_mod_cubic
    }

    /// Returns `gcd(f(x), x^q - x)` computed via the quotient-ring reduction.
    pub fn gcd(&self) -> &DensePolynomial<F> {
        &self.gcd
    }

    /// Returns whether `E(F_q)` has rational `2`-torsion.
    pub fn has_rational_two_torsion(&self) -> bool {
        self.has_rational_two_torsion
    }

    /// Returns `tr(π_q) mod 2`.
    pub fn trace_mod_2(&self) -> u8 {
        self.trace_mod_2
    }

    /// Returns the congruence `t ≡ trace_mod_2 (mod 2)`.
    pub fn trace_congruence(&self) -> Congruence {
        Congruence::new(BigUint::from(self.trace_mod_2), BigUint::from(2u8))
            .expect("mod 2 should define a valid congruence")
    }

    /// Returns whether `#E(F_q)` is even.
    ///
    /// Over odd `q`, the parity of `#E(F_q)` agrees with `tr(π_q) mod 2`
    /// because `q + 1` is even.
    pub fn group_order_is_even(&self) -> bool {
        self.has_rational_two_torsion
    }
}
