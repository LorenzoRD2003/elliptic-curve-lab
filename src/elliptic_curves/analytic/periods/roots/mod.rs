//! Validated cubic-root triples for the analytic Weierstrass model.
//!
//! This submodule is the local owner of the value
//! `WeierstrassCubicRoots` and of the first geometric questions one can ask
//! once those roots have been recovered:
//!
//! - are they well separated or nearly repeated?
//! - do they look approximately real, approximately conjugate, or fully
//!   generic complex?
//! - which Legendre reduction do they induce?
//!
//! Accordingly:
//!
//! - `value.rs` stores the validated triple and its symmetric-function data.
//! - `classification.rs` records the geometric/separation reports attached to
//!   one triple.
//!
//! This keeps the “shape of the roots” story with the roots themselves rather
//! than scattering it across the broader `periods` namespace.
mod classification;
mod value;

pub use classification::{
    CubicRootConfiguration, CubicRootConfigurationReport, CubicRootSeparation,
};
pub use value::WeierstrassCubicRoots;
