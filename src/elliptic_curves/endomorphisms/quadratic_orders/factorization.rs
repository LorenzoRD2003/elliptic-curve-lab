use num_bigint::{BigInt, BigUint};
use num_traits::{One, Signed};

use crate::elliptic_curves::{
    endomorphisms::quadratic_orders::{
        ImaginaryQuadraticOrder, ImaginaryQuadraticOrderError, QuadraticDiscriminant,
        QuadraticDiscriminantFactorizationError,
    },
    frobenius::FrobeniusDiscriminant,
};
use crate::numerics::quadratic_radicands::{QuadraticIntegerBasisKind, split_square_part};

/// Canonical factorization of a quadratic discriminant `Δ = v^2 D_K`,
/// where `D_K` is a fundamental discriminant.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuadraticDiscriminantFactorization {
    discriminant: QuadraticDiscriminant,
    conductor: BigUint,
    fundamental_discriminant: QuadraticDiscriminant,
}

impl QuadraticDiscriminant {
    /// Factors a negative quadratic-order discriminant as `Δ = v^2 D_K`,
    /// where `D_K` is the fundamental discriminant of an imaginary quadratic
    /// field and `v >= 1` is the positive square root factor.
    ///
    /// Aaccepts only negative discriminants congruent to `0` or `1` modulo `4`.
    ///
    /// Complexity: Dominated by `num-prime`.
    pub fn factorization(
        &self,
    ) -> Result<QuadraticDiscriminantFactorization, QuadraticDiscriminantFactorizationError> {
        self.validate_imaginary_discriminant()?;

        if self.is_congruent_to_1_mod_4() {
            self.factorization_for_one_mod_four()
        } else if self.is_congruent_to_0_mod_4() {
            self.factorization_for_zero_mod_four()
        } else {
            Err(QuadraticDiscriminantFactorizationError::InvalidQuadraticOrderDiscriminant)
        }
    }

    fn validate_imaginary_discriminant(
        &self,
    ) -> Result<(), QuadraticDiscriminantFactorizationError> {
        if self.is_zero() {
            return Err(QuadraticDiscriminantFactorizationError::ZeroDiscriminant);
        }
        if self.is_positive() {
            return Err(QuadraticDiscriminantFactorizationError::PositiveDiscriminant);
        }
        Ok(())
    }

    fn factorization_for_one_mod_four(
        &self,
    ) -> Result<QuadraticDiscriminantFactorization, QuadraticDiscriminantFactorizationError> {
        let (square_root_factor, squarefree_part) =
            split_square_part(&absolute_biguint(self.value()));
        let fundamental_discriminant = QuadraticDiscriminant::new(-BigInt::from(squarefree_part));

        Ok(self.build_factorization(square_root_factor, fundamental_discriminant))
    }

    fn factorization_for_zero_mod_four(
        &self,
    ) -> Result<QuadraticDiscriminantFactorization, QuadraticDiscriminantFactorizationError> {
        let quarter = self.value() / 4u8;
        let (square_root_factor, squarefree_part) = split_square_part(&absolute_biguint(&quarter));
        let signed_squarefree_part = -BigInt::from(squarefree_part);
        let basis_kind = QuadraticIntegerBasisKind::for_squarefree_part(&signed_squarefree_part);
        let conductor = basis_kind.discriminant_square_root_factor(square_root_factor);
        let fundamental_discriminant = QuadraticDiscriminant::new(
            basis_kind.fundamental_discriminant_for(&signed_squarefree_part),
        );

        Ok(self.build_factorization(conductor, fundamental_discriminant))
    }

    fn build_factorization(
        &self,
        conductor: BigUint,
        fundamental_discriminant: QuadraticDiscriminant,
    ) -> QuadraticDiscriminantFactorization {
        debug_assert!(fundamental_discriminant.is_fundamental());
        QuadraticDiscriminantFactorization {
            discriminant: self.clone(),
            conductor,
            fundamental_discriminant,
        }
    }
}

impl QuadraticDiscriminantFactorization {
    /// Factors one Frobenius discriminant `Δ_π = t^2 - 4q` as `Δ_π = v^2 D_K`.
    ///
    /// Complexity: dominated by `num-prime`.
    pub fn from_frobenius_discriminant(
        frobenius_discriminant: &FrobeniusDiscriminant,
    ) -> Result<Self, QuadraticDiscriminantFactorizationError> {
        frobenius_discriminant
            .quadratic_discriminant()
            .factorization()
    }

    /// Returns the original discriminant `Δ`.
    pub fn discriminant(&self) -> &QuadraticDiscriminant {
        &self.discriminant
    }

    /// Returns the positive square root factor `v`.
    ///
    /// For quadratic-order discriminants this is the conductor in the
    /// decomposition `Δ = v^2 D_K`.
    pub fn conductor(&self) -> &BigUint {
        &self.conductor
    }

    /// Returns the fundamental discriminant `D_K`.
    pub fn fundamental_discriminant(&self) -> &QuadraticDiscriminant {
        &self.fundamental_discriminant
    }

    /// Returns whether the factorization is trivial, i.e. `v = 1`.
    pub fn is_fundamental_already(&self) -> bool {
        self.conductor == BigUint::one()
    }

    /// Returns the maximal order `O_K`, corresponding to conductor `f = 1`.
    ///
    /// If `Δ = v^2 D_K`, then this helper discards the conductor `v` and
    /// returns the maximal order with discriminant `D_K`.
    ///
    /// Complexity: `Θ(1)`
    pub fn maximal_order(&self) -> Result<ImaginaryQuadraticOrder, ImaginaryQuadraticOrderError> {
        ImaginaryQuadraticOrder::new(self.fundamental_discriminant.clone(), BigUint::one())
    }
}

fn absolute_biguint(value: &BigInt) -> BigUint {
    value
        .abs()
        .to_biguint()
        .expect("absolute value of a non-zero integer should be non-negative")
}
