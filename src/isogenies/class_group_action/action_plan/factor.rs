use num_bigint::BigInt;
use num_traits::Zero;

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::BinaryQuadraticForm, quadratic_ideals::PrimeNormIdeal,
};

/// One nonzero local factor in an algebraic class-group action plan.
///
/// The factor records the prime-norm ideal selected by the caller, the reduced
/// form class associated to that ideal by the current ideal-to-form convention,
/// and the signed exponent found by the finite subgroup search.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ClassGroupActionPlanFactor {
    ideal: PrimeNormIdeal,
    generator_form: BinaryQuadraticForm,
    exponent: BigInt,
}

impl ClassGroupActionPlanFactor {
    pub(super) fn from_nonzero_exponents(
        local_ideals: &[PrimeNormIdeal],
        generator_forms: Vec<BinaryQuadraticForm>,
        exponents: Vec<BigInt>,
    ) -> Vec<Self> {
        local_ideals
            .iter()
            .cloned()
            .zip(generator_forms)
            .zip(exponents)
            .filter_map(|((ideal, generator_form), exponent)| {
                (!exponent.is_zero()).then(|| Self {
                    ideal,
                    generator_form,
                    exponent,
                })
            })
            .collect()
    }

    /// Returns the local prime-norm ideal used by this factor.
    pub(crate) fn ideal(&self) -> &PrimeNormIdeal {
        &self.ideal
    }

    /// Returns the reduced form class associated to the local ideal.
    pub(crate) fn generator_form(&self) -> &BinaryQuadraticForm {
        &self.generator_form
    }

    /// Returns the signed exponent assigned to this local generator.
    pub(crate) fn exponent(&self) -> &BigInt {
        &self.exponent
    }
}
