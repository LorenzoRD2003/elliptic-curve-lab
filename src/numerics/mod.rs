//! Shared numerical helpers for approximate and exact workflows.

mod bernoulli;
pub mod chinese_remainder;
mod comparison;
mod complex_algebra;
mod complex_path;
pub mod cornacchia;
mod euclidean_division;
mod gcd;
pub(crate) mod hensel;
pub(crate) mod integer_arithmetic;
mod linear_recurrence;
mod number_theory;
#[cfg(test)]
pub(crate) mod perfect_powers;
mod prime_powers;
pub mod quadratic_forms;
pub(crate) mod quadratic_radicands;
mod rational_arithmetic;
mod sigma;
mod simpson;
mod tolerance;

pub use bernoulli::bernoulli_number;
pub(crate) use comparison::{
    ComplexApproxComparison, ComplexDifferenceReport, HasComplexApproxComparison,
};
pub(crate) use complex_algebra::{
    cube_root_branches, is_near_pure_cubic_regime, primitive_cube_root_of_unity,
};
pub use complex_path::{ComplexLineSegment, ComplexRay};
pub(crate) use euclidean_division::{ceil_div_bigint_by_positive, floor_div_bigint_by_positive};
pub(crate) use gcd::{extended_gcd_bigint, gcd_bigint, gcd_biguint, inverse_mod_biguint};
pub(crate) use integer_arithmetic::{
    lcm_bigint, lcm_biguint, lcm_usize, pow_bigint_usize, quotients_by_distinct_prime_factors,
};
pub use linear_recurrence::OrderTwoLinearRecurrence;
pub use number_theory::{PositivePrimeError, is_squarefree, positive_divisors};
pub(crate) use number_theory::{
    distinct_prime_factors, exact_square_root, validate_positive_prime, valuation_biguint,
};
pub(crate) use prime_powers::{NormalizedPrimePowerFactorization, PrimePowerTable};
pub(crate) use rational_arithmetic::rational_denominator_power_clearance;
pub use sigma::{sigma_power_sum_factorized, sigma_power_sums_up_to};
pub use simpson::{
    SimpsonIntegrationError, SimpsonQuadratureDomain, SimpsonQuadratureDomainError,
    composite_simpson_integrate_complex_in_domain,
    composite_simpson_integrate_complex_simple_in_domain,
};
pub use tolerance::ApproxTolerance;
pub(crate) use tolerance::{
    projective_unit_singularity_distance, reciprocal_singularity_threshold,
};
