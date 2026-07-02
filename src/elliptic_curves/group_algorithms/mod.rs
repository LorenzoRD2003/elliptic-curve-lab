//! Curve-group algorithms that are not specific to one Weierstrass model.

mod cyclic_primary_order;
pub mod cyclic_roots;
mod small_finite_api;

pub(crate) use cyclic_primary_order::CyclicPrimaryOrderGroupCurveModel;
pub(crate) use small_finite_api::{
    group_exponent_by as shared_group_exponent_by, point_order_by as shared_point_order_by,
};
