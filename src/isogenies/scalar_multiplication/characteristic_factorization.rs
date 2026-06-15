use crate::numerics::{pow_u64_as_usize, square_u64_as_usize};

/// Characteristic-side factorization data for the scalar `n` in `[n]`.
///
/// The current public interpretation writes `n = p^e m`,
/// where `p` is the base-field characteristic and `gcd(m, p) = 1`.
///
/// This package then records:
///
/// - the exponent `e`
/// - the prime-to-characteristic factor `m`
/// - the degree `m^2` of the visible reduced factor currently exposed
/// - the residual characteristic-`p` degree bucket `p^(2e)`
///
/// Scope note:
/// the current crate does not yet refine that residual `p`-power contribution
/// into the finer ordinary/supersingular geometric subcases.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScalarCharacteristicFactorization {
    p_power_exponent: u32,
    separable_part: u64,
    separable_degree: usize,
    infinitesimal_degree: usize,
}

impl ScalarCharacteristicFactorization {
    /// Returns the exponent `e` in the factorization `n = p^e m`
    /// where `p` is the base-field characteristic and `gcd(m, p) = 1`.
    pub fn p_power_exponent(&self) -> u32 {
        self.p_power_exponent
    }

    /// Returns the prime-to-characteristic factor `m` in
    /// `n = p^e m`.
    pub fn separable_part(&self) -> u64 {
        self.separable_part
    }

    /// Returns the degree currently attributed to the visible reduced part.
    ///
    /// In the current implementation this is `m^2`, where `m = separable_part()`.
    pub fn separable_degree(&self) -> usize {
        self.separable_degree
    }

    /// Returns the residual characteristic-`p` degree bucket.
    ///
    /// In the current implementation this is `p^(2e)`, where `n = p^e m`.
    ///
    /// Scope note:
    /// this is the degree of the full `p`-power contribution, which the
    /// current public kernel description does not yet refine into ordinary
    /// versus supersingular subcases.
    pub fn infinitesimal_degree(&self) -> usize {
        self.infinitesimal_degree
    }

    pub(super) fn new(
        p_power_exponent: u32,
        separable_part: u64,
        separable_degree: usize,
        infinitesimal_degree: usize,
    ) -> Self {
        Self {
            p_power_exponent,
            separable_part,
            separable_degree,
            infinitesimal_degree,
        }
    }

    pub(super) fn from_scalar_and_characteristic(scalar: u64, characteristic: u64) -> Self {
        if characteristic == 0 || characteristic == 1 {
            return Self::new(0, scalar, square_u64_as_usize(scalar), 1);
        }

        let mut exponent = 0u32;
        let mut separable_part = scalar;
        while separable_part.is_multiple_of(characteristic) {
            separable_part /= characteristic;
            exponent += 1;
        }

        let separable_degree = square_u64_as_usize(separable_part);
        let infinitesimal_degree = pow_u64_as_usize(characteristic, exponent.saturating_mul(2));

        Self::new(
            exponent,
            separable_part,
            separable_degree,
            infinitesimal_degree,
        )
    }
}
