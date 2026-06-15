//! Public kernel-description surface for isogenies.
//!
//! This submodule separates four small stories:
//!
//! - `kernel_description.rs`: the top-level public enum and shared queries
//! - `reduced.rs`: fully reduced point-visible kernel data
//! - `nonreduced.rs`: purely infinitesimal kernel data
//! - `mixed.rs`: kernels with both reduced and infinitesimal parts

mod kernel_description;
mod mixed;
mod nonreduced;
mod reduced;

pub use kernel_description::KernelDescription;
pub use mixed::MixedKernelDescription;
pub use nonreduced::NonReducedKernelDescription;
pub use reduced::ReducedKernelDescription;
