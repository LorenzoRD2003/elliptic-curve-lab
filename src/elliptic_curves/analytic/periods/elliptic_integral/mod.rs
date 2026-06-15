//! Complete elliptic integrals and AGM-side period data.
//!
//! After choosing a Legendre parameter `λ`, the normalized half-period story
//! is governed by the complete elliptic integrals `K(λ)` and `K(1-λ)`.
//! In the current crate those values are approximated through the complex
//! arithmetic-geometric mean.
//!
//! The submodule split follows that narrative:
//!
//! - `agm.rs` contains the raw complex AGM primitive and its trace/report
//!   types.
//! - `k.rs` packages the resulting approximations as
//!   complete-elliptic-integral reports.
//! - `period_report.rs` combines the two `K` values into the Legendre-side
//!   period report and resulting `τ` candidate.
//!
//! The current APIs are still educational: they expose branch choices,
//! iteration counts, and numerical status rather than hiding the AGM route
//! behind one opaque complex number.
mod agm;
mod k;
mod period_report;

pub use agm::{
    ComplexAgmBranchChoice, ComplexAgmConfig, ComplexAgmIteration, ComplexAgmResult,
    ComplexAgmStatus, ComplexAgmTrace, complex_agm, complex_agm_trace,
};
pub use k::{CompleteEllipticIntegralKApprox, CompleteEllipticIntegralKMetadata};
pub use period_report::LegendrePeriodIntegralReport;
