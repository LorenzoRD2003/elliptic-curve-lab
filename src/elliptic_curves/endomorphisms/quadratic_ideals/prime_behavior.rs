use num_bigint::BigUint;
use num_traits::Zero;

use crate::elliptic_curves::endomorphisms::{
    quadratic_ideals::QuadraticPrimeBehaviorError, quadratic_orders::ImaginaryQuadraticOrder,
};
use crate::numerics::{
    hensel::sqrt_mod_odd_prime_power, positive_mod_biguint, validate_positive_prime,
};

/// Local behavior of a prime `ℓ` in an imaginary quadratic order `O_f`.
///
/// This value is not yet an ideal, an ideal class, or an action on curves. It
/// is the first local arithmetic classification needed before those later
/// objects can be introduced honestly.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum QuadraticPrimeBehavior {
    /// `ℓ` splits and gives two local square roots of `Δ` modulo `ℓ`.
    ///
    /// The tuple stores canonical representatives `(r₁, r₂)` with `r₁ < r₂`.
    /// They satisfy `rᵢ² ≡ Δ (mod ℓ)` and `r₁ + r₂ = ℓ`.
    ///
    /// In the later volcano story, this is the case compatible with two
    /// horizontal directions of degree `ℓ`, once the graph-side evidence has
    /// also certified horizontality.
    Split { roots: (BigUint, BigUint) },
    /// `ℓ` is inert, equivalently `(Δ/ℓ) = -1`.
    ///
    /// This is the local obstruction to horizontal `ℓ`-isogeny motion.
    Inert,
    /// `ℓ` ramifies, equivalently `(Δ/ℓ) = 0` in the invertible-prime regime.
    Ramified { root: BigUint },
    /// `ℓ | f`, so `ℓ` is not invertible in the non-maximal order `O_f`.
    NonInvertibleBecauseDividesConductor,
}

impl ImaginaryQuadraticOrder {
    /// Classifies the behavior of an odd prime `ℓ` in `O_f`.
    ///
    /// The method first checks whether `ℓ | f`. In that case it returns
    /// [`QuadraticPrimeBehavior::NonInvertibleBecauseDividesConductor`], because
    /// the prime is not in the invertible ideal regime of the non-maximal order.
    ///
    /// For odd `ℓ ∤ f`, the order discriminant `Δ = f² D_K` is a unit modulo
    /// `ℓ` unless `ℓ` ramifies in the quadratic field. The returned cases
    /// correspond to the Kronecker-Legendre value `(Δ/ℓ)`:
    ///
    /// - `1`: split, with the two roots of `x² ≡ Δ (mod ℓ)`;
    /// - `0`: ramified, with the repeated root `0`;
    /// - `-1`: inert.
    ///
    /// The first slice deliberately rejects `ℓ = 2` unless `2 | f`: the
    /// dyadic splitting rules depend on `Δ mod 8`, while the current
    /// [`QuadraticPrimeBehavior::Split`] variant stores two roots modulo `ℓ`.
    ///
    /// Complexity: prime validation is dominated by `num-prime`; the split
    /// case is dominated by the crate's Tonelli-Shanks square-root route.
    pub fn prime_behavior(
        &self,
        ell: &BigUint,
    ) -> Result<QuadraticPrimeBehavior, QuadraticPrimeBehaviorError> {
        validate_positive_prime(ell)?;

        if (self.conductor() % ell).is_zero() {
            return Ok(QuadraticPrimeBehavior::NonInvertibleBecauseDividesConductor);
        } else if ell == &BigUint::from(2u8) {
            return Err(QuadraticPrimeBehaviorError::UnsupportedPrimeTwo);
        }

        let discriminant_mod_ell = positive_mod_biguint(self.discriminant().value(), ell);
        if discriminant_mod_ell.is_zero() {
            return Ok(QuadraticPrimeBehavior::Ramified {
                root: BigUint::zero(),
            });
        }

        match sqrt_mod_odd_prime_power(self.discriminant().value(), ell, 1) {
            Ok((left, right)) => Ok(QuadraticPrimeBehavior::Split {
                roots: (left, right),
            }),
            Err(_) => Ok(QuadraticPrimeBehavior::Inert),
        }
    }
}
