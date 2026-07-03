use crate::elliptic_curves::short_weierstrass::rational_torsion::RationalTorsionError;

#[test]
fn staged_error_variants_have_readable_messages() {
    assert!(
        RationalTorsionError::IntegralModelUnavailable
            .to_string()
            .contains("integral short-Weierstrass model")
    );
    assert!(
        RationalTorsionError::InconsistentMazurShape { point_count: 11 }
            .to_string()
            .contains("11 torsion points")
    );
}
