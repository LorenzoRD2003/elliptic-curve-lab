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
    assert!(
        RationalTorsionError::InconsistentReportGroup {
            group_cardinality: 6,
            point_count: 5,
        }
        .to_string()
        .contains("5 points")
    );
    assert!(
        RationalTorsionError::InvalidCandidateAccounting {
            candidate_count: 2,
            point_count: 3,
        }
        .to_string()
        .contains("2 candidates")
    );
}
