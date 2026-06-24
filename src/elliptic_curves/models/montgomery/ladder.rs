//! Montgomery ladder on the `x`-line of the normalized model.
//!
//! This module implements the first educational ladder surface for
//! `v^2 = x^3 + A x^2 + x`.
//!
//! The maintained state is one neighboring pair
//!
//! - `R0 = [k]P`
//! - `R1 = [k+1]P`
//!
//! so the fixed differential invariant is
//!
//! - `R1 - R0 = P`.
//!
//! The ladder starts from `(R0, R1) = (O, P)`, which already satisfies that
//! invariant. Each processed scalar bit updates the pair using only the Stage B
//! differential primitives:
//!
//! - if the next bit is `0`, the new pair is
//!   `([2k]P, [2k+1]P) = ([2]R0, R0 + R1)`
//! - if the next bit is `1`, the new pair is
//!   `([2k+1]P, [2k+2]P) = (R0 + R1, [2]R1)`
//!
//! In both branches, the new difference is still `P`, because
//!
//! - `[2k+1]P - [2k]P = P`
//! - `[2k+2]P - [2k+1]P = P`.
//!
//! Reading the scalar bits from most significant to least significant therefore
//! evolves `k` by the usual binary recurrence
//!
//! - `k -> 2k` when the bit is `0`
//! - `k -> 2k + 1` when the bit is `1`
//!
//! while never requiring the sign of `y(P)`.

use crate::elliptic_curves::{MontgomeryCurve, MontgomeryXzPoint, NormalizedMontgomeryCurve};
use crate::fields::traits::{Field, SqrtField};

impl<F: Field> NormalizedMontgomeryCurve<F> {
    fn ladder_state_from_base(
        &self,
        base: &MontgomeryXzPoint<F>,
        scalar: u64,
    ) -> (MontgomeryXzPoint<F>, MontgomeryXzPoint<F>) {
        let mut r0 = MontgomeryXzPoint::Infinity;
        let mut r1 = base.clone();

        for bit_index in (0..u64::BITS).rev() {
            if ((scalar >> bit_index) & 1) == 0 {
                let (next_r0, next_r1) = self
                    .x_dbl_add(&r0, &r1, base)
                    .expect("the ladder invariant keeps x(P) as the difference witness");
                r0 = next_r0;
                r1 = next_r1;
            } else {
                let (next_r1, next_r0) = self
                    .x_dbl_add(&r1, &r0, base)
                    .expect("the ladder invariant keeps x(P) as the difference witness");
                r0 = next_r0;
                r1 = next_r1;
            }
        }

        (r0, r1)
    }

    /// Returns the final ladder pair `(x([n]P), x([n+1]P))` from the affine
    /// `x`-coordinate of `P`.
    ///
    /// The maintained invariant is that the two projective `x`-line values
    /// differ by the fixed base point `P`, so each bit step can use only the
    /// differential primitives `xDBL` and `xADD`.
    ///
    /// Complexity: `Θ(log n)` differential steps for a `u64` scalar `n`.
    pub(crate) fn ladder_xz_pair(
        &self,
        base_x: F::Elem,
        scalar: u64,
    ) -> (MontgomeryXzPoint<F>, MontgomeryXzPoint<F>) {
        self.ladder_state_from_base(&MontgomeryXzPoint::from_affine_x(base_x), scalar)
    }

    /// Computes the `x`-line value of `[n]P` from the affine `x`-coordinate of
    /// `P`, without requiring the sign of `y(P)`.
    ///
    /// This is the first educational Montgomery ladder surface in the repo.
    /// It follows a fixed bit-by-bit schedule, but it is not claimed to be a
    /// production constant-time implementation across all backends.
    ///
    /// Complexity: `Θ(log n)` differential steps for a `u64` scalar `n`.
    pub fn ladder_x(&self, base_x: F::Elem, scalar: u64) -> MontgomeryXzPoint<F> {
        self.ladder_xz_pair(base_x, scalar).0
    }
}

impl<F: Field + SqrtField> MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    /// Computes the `x`-line value of `[n]P` from the affine `x`-coordinate of
    /// `P` when the same-field normalization to `B = 1` is available.
    ///
    /// Complexity: the same `Θ(log n)` ladder cost as the normalized route,
    /// plus one normalization lookup for the current curve descriptor.
    pub fn try_ladder_x(
        &self,
        x: F::Elem,
        scalar: u64,
    ) -> Result<MontgomeryXzPoint<F>, crate::elliptic_curves::MontgomeryNormalizationError> {
        self.try_as_normalized_montgomery()
            .map(|normalized| normalized.ladder_x(x, scalar))
    }
}
