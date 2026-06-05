//! Shared numerical helpers for approximate and exact workflows.

pub mod bernoulli;
pub mod sigma;
pub mod tolerance;

pub use bernoulli::bernoulli_number;
pub use sigma::{sigma_power_sum_factorized, sigma_power_sum_naive, sigma_power_sums_up_to};
pub use tolerance::ApproxTolerance;
