use super::*;

#[test]
fn config_constructor_preserves_caller_supplied_values() {
    let tolerance = ApproxTolerance::new(1.0e-8, 2.0e-8);
    let config = PeriodRecoveryConfig::new(tolerance, 9, 7, 192, 3, 11).unwrap();

    assert_eq!(config.tolerance(), tolerance);
    assert_eq!(config.newton_max_iterations(), 9);
    assert_eq!(config.agm_max_iterations(), 7);
    assert_eq!(config.abel_jacobi_integration_steps(), 192);
    assert_eq!(config.branch_lattice_search_radius(), 3);
    assert_eq!(config.fundamental_domain_reduction_max_steps(), 11);
}

#[test]
fn config_rejects_zero_budgets() {
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 0, 1, 1, 1, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 0, 1, 1, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 1, 0, 1, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 1, 1, 0, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 1, 1, 1, 0),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
}

#[test]
fn config_presets_are_ordered_and_explicit() {
    let educational = PeriodRecoveryConfig::educational_default();
    let strict = PeriodRecoveryConfig::strict();
    let loose = PeriodRecoveryConfig::loose();

    assert_eq!(
        educational.tolerance(),
        ApproxTolerance::educational_default()
    );
    assert_eq!(educational.newton_max_iterations(), 12);
    assert_eq!(educational.agm_max_iterations(), 10);
    assert_eq!(educational.abel_jacobi_integration_steps(), 256);
    assert_eq!(educational.branch_lattice_search_radius(), 2);
    assert_eq!(educational.fundamental_domain_reduction_max_steps(), 16);

    assert_eq!(strict.tolerance(), ApproxTolerance::strict());
    assert!(strict.newton_max_iterations() > educational.newton_max_iterations());
    assert!(strict.agm_max_iterations() > educational.agm_max_iterations());
    assert!(strict.abel_jacobi_integration_steps() > educational.abel_jacobi_integration_steps());
    assert!(strict.branch_lattice_search_radius() > educational.branch_lattice_search_radius());
    assert!(
        strict.fundamental_domain_reduction_max_steps()
            > educational.fundamental_domain_reduction_max_steps()
    );

    assert_eq!(loose.tolerance(), ApproxTolerance::loose());
    assert!(loose.newton_max_iterations() < educational.newton_max_iterations());
    assert!(loose.agm_max_iterations() < educational.agm_max_iterations());
    assert!(loose.abel_jacobi_integration_steps() < educational.abel_jacobi_integration_steps());
    assert!(loose.branch_lattice_search_radius() < educational.branch_lattice_search_radius());
    assert!(
        loose.fundamental_domain_reduction_max_steps()
            < educational.fundamental_domain_reduction_max_steps()
    );
}
