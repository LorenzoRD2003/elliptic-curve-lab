use crate::elliptic_curves::endomorphisms::{
    ImaginaryQuadraticOrderError, QuadraticDiscriminant, QuadraticDiscriminantFactorization,
    QuadraticOrderIndexError,
};
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Zero};

/// Imaginary quadratic order `O_f = ℤ + f O_K`.
///
/// Here `D_K < 0` is the fundamental discriminant of the imaginary quadratic
/// field `K`, `f >= 1` is the conductor, and the relation is `Δ(O_f) = f^2 D_K`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImaginaryQuadraticOrder {
    fundamental_discriminant: QuadraticDiscriminant,
    conductor: BigUint,
    discriminant: QuadraticDiscriminant,
}

impl ImaginaryQuadraticOrder {
    /// Builds the imaginary quadratic order `O_f = ℤ + f O_K` from `D_K` and `f`.
    ///
    /// Complexity: `Θ(1)` big-integer arithmetic.
    pub fn new(
        fundamental_discriminant: QuadraticDiscriminant,
        conductor: BigUint,
    ) -> Result<Self, ImaginaryQuadraticOrderError> {
        validate_imaginary_order_inputs(&fundamental_discriminant, &conductor)?;
        let conductor_bigint = BigInt::from_biguint(Sign::Plus, conductor.clone());
        let conductor_squared = &conductor_bigint * &conductor_bigint;
        let discriminant =
            QuadraticDiscriminant::new(conductor_squared * fundamental_discriminant.value());

        Ok(Self {
            fundamental_discriminant,
            conductor,
            discriminant,
        })
    }

    /// Builds the order from a canonical factorization `Δ = f^2 D_K`.
    ///
    /// Complexity: `Θ(1)` once the factorization has already been computed.
    pub fn from_factorization(
        factorization: QuadraticDiscriminantFactorization,
    ) -> Result<Self, ImaginaryQuadraticOrderError> {
        Self::new(
            factorization.fundamental_discriminant().clone(),
            factorization.conductor().clone(),
        )
    }

    /// Builds the order from a discriminant `Δ` by first factoring it as `Δ = f^2 D_K`.
    ///
    /// Complexity: dominated by `num-prime`.
    pub fn from_discriminant(
        discriminant: &QuadraticDiscriminant,
    ) -> Result<Self, ImaginaryQuadraticOrderError> {
        let factorization = discriminant
            .factorization()
            .map_err(|_| ImaginaryQuadraticOrderError::NonImaginaryOrderDiscriminant)?;
        Self::from_factorization(factorization)
    }

    /// Returns the fundamental discriminant `D_K`.
    pub fn fundamental_discriminant(&self) -> &QuadraticDiscriminant {
        &self.fundamental_discriminant
    }

    /// Returns the conductor `f`.
    pub fn conductor(&self) -> &BigUint {
        &self.conductor
    }

    /// Returns the order discriminant `Δ(O_f) = f^2 D_K`.
    pub fn discriminant(&self) -> &QuadraticDiscriminant {
        &self.discriminant
    }

    /// Returns whether the order is maximal, equivalently whether `f = 1`.
    pub fn is_maximal(&self) -> bool {
        self.conductor == BigUint::one()
    }

    /// Returns whether `self` and `other` lie in the same imaginary quadratic field.
    pub fn same_quadratic_field(&self, other: &Self) -> bool {
        self.fundamental_discriminant == other.fundamental_discriminant
    }

    /// Returns whether `self` is a suborder of `other`.
    ///
    /// For orders inside the same imaginary quadratic field, one has
    /// `O_f ⊆ O_g` if and only if `g | f`.
    pub fn is_suborder_of(&self, other: &Self) -> bool {
        self.same_quadratic_field(other) && (&self.conductor % &other.conductor).is_zero()
    }

    /// Returns whether `self` is an overorder of `other`.
    ///
    /// Equivalently, this is the reverse relation of [`Self::is_suborder_of`].
    pub fn is_overorder_of(&self, other: &Self) -> bool {
        other.is_suborder_of(self)
    }

    /// Returns the index `[self : suborder]` when `suborder ⊆ self`.
    ///
    /// For imaginary quadratic orders in the same field, if `self = O_{f_2}` and
    /// `suborder = O_{f_1}` with `f_2 | f_1`, then `[self : suborder] = f_1 / f_2`.
    ///
    /// Complexity: `Θ(1)` big-integer arithmetic.
    pub fn index_of_suborder(&self, suborder: &Self) -> Result<BigUint, QuadraticOrderIndexError> {
        if !self.same_quadratic_field(suborder) {
            return Err(QuadraticOrderIndexError::DifferentQuadraticFields);
        }
        if !suborder.is_suborder_of(self) {
            return Err(QuadraticOrderIndexError::NotSuborder);
        }

        Ok(&suborder.conductor / &self.conductor)
    }
}

fn validate_imaginary_order_inputs(
    fundamental_discriminant: &QuadraticDiscriminant,
    conductor: &BigUint,
) -> Result<(), ImaginaryQuadraticOrderError> {
    if !fundamental_discriminant.is_negative() {
        return Err(ImaginaryQuadraticOrderError::NonNegativeFundamentalDiscriminant);
    }
    if !fundamental_discriminant.is_fundamental() {
        return Err(ImaginaryQuadraticOrderError::NonFundamentalDiscriminant);
    }
    if conductor.is_zero() {
        return Err(ImaginaryQuadraticOrderError::ZeroConductor);
    }
    Ok(())
}
