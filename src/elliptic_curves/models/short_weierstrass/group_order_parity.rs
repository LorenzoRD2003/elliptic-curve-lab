use std::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    frobenius::{
        HasseInterval,
        hasse::search::{HasseBsgsConfig, HasseBsgsParity},
    },
    traits::HasseIntervalSearchCurveModel,
};
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

/// Parity of the finite group order `#E(F_q)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum GroupOrderParity {
    Even,
    Odd,
}

impl GroupOrderParity {
    #[allow(dead_code)]
    pub(crate) fn is_even(self) -> bool {
        matches!(self, Self::Even)
    }
}

/// Returns the parity of `#E(F_q)` by testing for rational `2`-torsion.
///
/// For a short-Weierstrass curve `E : y^2 = x^3 + ax + b` over `F_q`, the
/// group order is even if and only if the cubic `f(x) = x^3 + ax + b` has
/// a root in `F_q`, equivalently if and only if
///
/// `deg gcd(x^q - x, f(x)) > 0`.
///
/// The current implementation computes this gcd criterion through the quotient
/// ring `F_q[x]/(f(x))`: it first computes `x^q mod f(x)` by repeated
/// squaring, then takes the gcd of `f(x)` with `(x^q mod f(x)) - x`.
/// This avoids materializing the degree-`q` polynomial `x^q - x`.
///
/// Complexity: `Θ(log q)` polynomial squarings/multiplications and reductions
/// in the quotient by a cubic, plus one gcd of degree at most `3`. Counting field
/// operations under the current dense backend, this is `Θ(log q)` field
/// operations with small constant factors.
impl<F: FiniteField> ShortWeierstrassCurve<F> {
    #[allow(dead_code)]
    pub(crate) fn group_order_parity_from_two_torsion(&self) -> GroupOrderParity {
        let cubic = self.to_cubic();
        let x = DensePolynomial::new(vec![F::zero(), F::one()]);
        let q = F::cardinality().unwrap();
        let x_q_mod_cubic = DensePolynomial::pow_mod(&x, q, &cubic)
            .expect("short-Weierstrass cubic is a non-zero modulus");
        let gcd = cubic.gcd(&x_q_mod_cubic.sub(&x));

        if gcd.degree().is_some_and(|degree| degree > 0) {
            GroupOrderParity::Even
        } else {
            GroupOrderParity::Odd
        }
    }

    /// Searches the discrete Hasse interval `H(q)` with BSGS after first
    /// determining the parity of `#E(F_q)` from rational `2`-torsion.
    ///
    /// This is the current public BSGS route for finding one annihilating
    /// multiple `M ∈ H(q)` with `[M]P = O` on the short-Weierstrass family.
    /// It is specific to the current model because it also uses the
    /// short-Weierstrass-only parity test from rational `2`-torsion.
    ///
    /// Complexity: one parity test of
    /// [`Self::group_order_parity_from_two_torsion`] plus a parity-restricted
    /// BSGS search on roughly half of `H(q)`. The resulting group-operation
    /// count is still `Θ(∜q)`, but with the baby-step and giant-step sizes
    /// both reduced by a factor of about `√2`.
    pub fn find_annihilating_multiple_in_hasse_interval_bsgs(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<Option<u128>, CurveError>
    where
        AffinePoint<F>: Clone + Eq + Hash,
    {
        let interval = HasseInterval::for_q(F::order())?;
        self.find_annihilating_multiple_in_interval_bsgs(point, interval)
    }

    /// Searches one Hasse interval with BSGS after first determining the
    /// parity of `#E(F_q)` from rational `2`-torsion.
    ///
    /// This is the curve-specific specialization of the generic Hasse BSGS
    /// engine: after deciding whether `#E(F_q)` is even or odd, it restricts
    /// the search to that parity class. Internally this amounts to advancing
    /// by `[2]P` in the baby-step side, matching the optimization described in
    /// the notes.
    ///
    /// Complexity: one parity test of
    /// [`Self::group_order_parity_from_two_torsion`] plus a parity-restricted
    /// BSGS search on roughly half of `H(q)`. The resulting group-operation
    /// count is still `Θ(∜q)`, but with the baby-step and giant-step sizes
    /// both reduced by a factor of about `√2`.
    pub(crate) fn find_annihilating_multiple_in_interval_bsgs(
        &self,
        point: &AffinePoint<F>,
        interval: HasseInterval,
    ) -> Result<Option<u128>, CurveError>
    where
        AffinePoint<F>: Clone + Eq + Hash,
    {
        let known_parity = match self.group_order_parity_from_two_torsion() {
            GroupOrderParity::Even => HasseBsgsParity::Even,
            GroupOrderParity::Odd => HasseBsgsParity::Odd,
        };
        HasseIntervalSearchCurveModel::find_annihilating_multiple_in_interval_bsgs_with_config(
            self,
            point,
            interval,
            HasseBsgsConfig::new().with_known_parity(known_parity),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::GroupOrderParity;
    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        traits::{EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel},
    };
    use crate::fields::{Fp, traits::Field};

    type F7 = Fp<7>;
    type F19 = Fp<19>;

    crate::fields::extension_field::define_fp_quadratic_extension!(
        spec: F19Sqrt2ParitySpec,
        field: F19Sqrt2Parity,
        base: F19,
        non_residue: 2,
        name: "F19(sqrt(2)) for group-order parity tests",
    );

    fn parity_from_order(order: usize) -> GroupOrderParity {
        if order.is_multiple_of(2) {
            GroupOrderParity::Even
        } else {
            GroupOrderParity::Odd
        }
    }

    #[test]
    fn two_torsion_parity_matches_small_prime_field_group_order() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");

        assert_eq!(
            curve.group_order_parity_from_two_torsion(),
            parity_from_order(curve.order())
        );
    }

    #[test]
    fn two_torsion_parity_detects_even_group_order_when_cubic_has_root() {
        let curve = ShortWeierstrassCurve::<F19>::new(F19::from_i64(-1), F19::zero())
            .expect("valid curve with x(x^2-1) cubic");

        assert_eq!(
            curve.group_order_parity_from_two_torsion(),
            GroupOrderParity::Even
        );
        assert_eq!(curve.order() % 2, 0);
    }

    #[test]
    fn two_torsion_parity_matches_small_extension_field_group_order() {
        let curve = ShortWeierstrassCurve::<F19Sqrt2Parity>::new(
            F19Sqrt2Parity::from_i64(2),
            F19Sqrt2Parity::from_i64(3),
        )
        .expect("valid extension-field curve");

        assert_eq!(
            curve.group_order_parity_from_two_torsion(),
            parity_from_order(curve.order())
        );
    }

    #[test]
    fn public_hasse_bsgs_helper_finds_an_annihilating_multiple() {
        let curve = ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("valid curve");
        let point = curve
            .generator()
            .expect("the sample curve should be cyclic");

        let multiple = curve
            .find_annihilating_multiple_in_hasse_interval_bsgs(&point)
            .expect("public Hasse BSGS helper should succeed")
            .expect("Hasse's theorem should guarantee an annihilating multiple");

        assert!(curve.is_torsion_point(
            &point,
            u64::try_from(multiple).expect("small Hasse multiple should fit in u64")
        ));
    }
}
