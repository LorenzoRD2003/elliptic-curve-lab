use crate::elliptic_curves::endomorphisms::{
    ImaginaryQuadraticOrder, ImaginaryQuadraticOrderError, QuadraticDiscriminant,
    QuadraticDiscriminantFactorizationError, QuadraticDiscriminantMod4,
};
use num_bigint::{BigInt, BigUint};
use num_prime::nt_funcs::factorize;
use num_traits::{One, Signed, Zero};

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

        match self.mod_4_class() {
            QuadraticDiscriminantMod4::One => self.factorization_for_one_mod_four(),
            QuadraticDiscriminantMod4::Zero => self.factorization_for_zero_mod_four(),
            QuadraticDiscriminantMod4::Other(_) => {
                Err(QuadraticDiscriminantFactorizationError::InvalidQuadraticOrderDiscriminant)
            }
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

        let (conductor, fundamental_discriminant) = factorization_data_from_even_squarefree_part(
            square_root_factor,
            signed_squarefree_part,
        )?;

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

fn split_square_part(value: &BigUint) -> (BigUint, BigUint) {
    if value.is_zero() {
        return (BigUint::one(), BigUint::zero());
    }

    let factors = factorize(value.clone());
    let mut square_root_factor = BigUint::one();
    let mut squarefree_part = BigUint::one();

    for (prime, exponent) in factors {
        let half = exponent / 2;
        if half > 0 {
            square_root_factor *= prime.pow(half as u32);
        }
        if exponent % 2 == 1 {
            squarefree_part *= prime;
        }
    }

    (square_root_factor, squarefree_part)
}

fn factorization_data_from_even_squarefree_part(
    square_root_factor: BigUint,
    signed_squarefree_part: BigInt,
) -> Result<(BigUint, QuadraticDiscriminant), QuadraticDiscriminantFactorizationError> {
    let squarefree_mod_four = normalized_mod_four(&signed_squarefree_part);

    if squarefree_mod_four == BigInt::one() {
        Ok((
            square_root_factor * BigUint::from(2u8),
            QuadraticDiscriminant::new(signed_squarefree_part),
        ))
    } else if squarefree_mod_four == BigInt::from(2u8) || squarefree_mod_four == BigInt::from(3u8) {
        Ok((
            square_root_factor,
            QuadraticDiscriminant::new(BigInt::from(4u8) * signed_squarefree_part),
        ))
    } else {
        Err(QuadraticDiscriminantFactorizationError::InvalidQuadraticOrderDiscriminant)
    }
}

fn absolute_biguint(value: &BigInt) -> BigUint {
    value
        .abs()
        .to_biguint()
        .expect("absolute value of a non-zero integer should be non-negative")
}

fn normalized_mod_four(value: &BigInt) -> BigInt {
    ((value % 4u8) + BigInt::from(4u8)) % 4u8
}
