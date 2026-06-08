//! Shared numerical helpers for approximate and exact workflows.

mod bernoulli;
mod comparison;
mod complex_path;
mod linear_recurrence;
mod number_theory;
mod sigma;
mod simpson;
mod tolerance;

pub use bernoulli::bernoulli_number;
pub(crate) use comparison::{
    ComplexApproxComparison, ComplexDifferenceReport, HasComplexApproxComparison,
};
pub use complex_path::{ComplexLineSegment, ComplexRay};
pub use linear_recurrence::OrderTwoLinearRecurrence;
pub use number_theory::{PositivePrimeError, is_squarefree, positive_divisors, valuation_biguint};
pub use sigma::{sigma_power_sum_factorized, sigma_power_sum_naive, sigma_power_sums_up_to};
pub use simpson::{
    SimpsonIntegrationError, SimpsonQuadratureDomain, SimpsonQuadratureDomainError,
    composite_simpson_integrate_complex_in_domain,
    composite_simpson_integrate_complex_simple_in_domain,
};
pub use tolerance::ApproxTolerance;
