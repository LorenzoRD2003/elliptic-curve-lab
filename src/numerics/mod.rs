//! Shared numerical helpers for approximate and exact workflows.

pub mod bernoulli;
pub mod complex_path;
pub mod sigma;
pub mod simpson;
pub mod tolerance;

pub use bernoulli::bernoulli_number;
pub use complex_path::{ComplexLineSegment, ComplexRay};
pub use sigma::{sigma_power_sum_factorized, sigma_power_sum_naive, sigma_power_sums_up_to};
pub use simpson::{
    SimpsonIntegrationError, SimpsonQuadratureDomain, SimpsonQuadratureDomainError,
    composite_simpson_integrate_complex_in_domain,
    composite_simpson_integrate_complex_simple_in_domain,
};
pub use tolerance::ApproxTolerance;
