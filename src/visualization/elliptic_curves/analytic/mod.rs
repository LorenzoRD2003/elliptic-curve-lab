mod curve;
mod formatting;
mod inverse_uniformization;
mod lattice;
mod modular;
mod periods;
#[cfg(test)]
mod tests;
mod torsion;

#[cfg(test)]
pub(crate) use formatting::format_complex_scalar_compact;

pub use curve::{
    describe_analytic_curve_membership, describe_analytic_invariants, describe_torus_to_curve_map,
    describe_weierstrass_differential_equation, describe_weierstrass_p_approx,
    describe_weierstrass_p_derivative_approx, format_analytic_cubic_model,
    format_short_weierstrass_over_complex,
};
pub use inverse_uniformization::{
    describe_invariant_recovery_validation_report,
    describe_inverse_uniformization_j_validation_report,
    describe_point_roundtrip_validation_config, describe_point_roundtrip_validation_report,
};
pub use lattice::{
    describe_complex_lattice, describe_eisenstein_sum, describe_q_parameter,
    describe_truncation_convergence,
};
pub use modular::{
    describe_fundamental_domain_reduction_report, describe_fundamental_domain_reduction_step,
    describe_j_invariant_comparison, describe_modular_invariance_report, describe_modular_matrix,
};
pub use periods::{
    describe_canonical_tau_recovery_report, describe_cubic_root_configuration_report,
    describe_cubic_root_recovery_report, describe_legendre_parameter,
    describe_legendre_parameter_conditioning, describe_legendre_parameter_orbit,
    describe_legendre_reduction, describe_legendre_reduction_report,
    describe_numerical_recovery_metadata, describe_period_basis_recovery_report,
    describe_period_lattice, describe_period_recovery_config, describe_period_recovery_report,
    describe_recovered_period_basis, describe_recovered_period_basis_report,
    describe_tau_recovery_report, describe_weierstrass_cubic_roots,
};
pub use torsion::{
    describe_analytic_division_polynomial_comparison,
    describe_analytic_even_division_polynomial_report,
    describe_analytic_odd_division_polynomial_report, describe_analytic_torsion_point_approx,
};
