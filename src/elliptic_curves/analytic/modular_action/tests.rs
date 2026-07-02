use num_bigint::BigInt;
use num_complex::Complex64;
use num_traits::One;
use proptest::prelude::*;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ApproxTolerance, LatticeSumTruncation, UpperHalfPlanePoint,
    modular_action::ModularMatrix,
};
use crate::proptest_support::{
    config::AnalyticStrategyConfig, elliptic_curves::arb_upper_half_plane_point,
};

fn modular_matrix_strategy() -> impl Strategy<Value = ModularMatrix> {
    ((-12i32..=12), (-12i32..=12), (-8i32..=8)).prop_filter_map(
        "pseudo-random coprime lower rows should produce a matrix in SL_2(ℤ)",
        |(c, d, shift)| {
            if c == 0 && d == 0 {
                return None;
            }

            let (gcd, bezout_c, bezout_d) = extended_gcd(c, d);
            if gcd != 1 {
                return None;
            }

            let a = bezout_d + shift * c;
            let b = -bezout_c + shift * d;

            ModularMatrix::new(a, b, c, d).ok()
        },
    )
}

fn modular_matrix_word_strategy() -> impl Strategy<Value = ModularMatrix> {
    prop::collection::vec(
        prop_oneof![Just(ModularMatrix::s()), Just(ModularMatrix::t())],
        0..8,
    )
    .prop_map(|word| {
        let mut product = ModularMatrix::identity();
        for generator in word {
            product = product
                .compose(&generator)
                .expect("products of modular generators stay in SL_2(ℤ)");
        }

        product
    })
}

fn modular_matrix_pair_strategy() -> impl Strategy<Value = (ModularMatrix, ModularMatrix)> {
    (modular_matrix_strategy(), modular_matrix_word_strategy())
        .prop_map(|(gamma, delta)| (gamma, delta))
}

fn extended_gcd(a: i32, b: i32) -> (i32, i32, i32) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (1, 0);
    let (mut old_t, mut t) = (0, 1);

    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
        (old_t, t) = (t, old_t - quotient * t);
    }

    if old_r < 0 {
        (-old_r, -old_s, -old_t)
    } else {
        (old_r, old_s, old_t)
    }
}

#[test]
fn constructor_accepts_standard_generators() {
    assert_eq!(ModularMatrix::identity().determinant(), BigInt::one());
    assert_eq!(ModularMatrix::s().determinant(), BigInt::one());
    assert_eq!(ModularMatrix::t().determinant(), BigInt::one());
}

#[test]
fn constructor_rejects_non_unimodular_matrix() {
    assert_eq!(
        ModularMatrix::new(1, 2, 3, 4),
        Err(AnalyticCurveError::InvalidModularMatrix)
    );
}

#[test]
fn accessors_match_constructor_entries() {
    let matrix = ModularMatrix::new(2, 1, 3, 2).expect("matrix should be valid");

    assert_eq!(matrix.a(), &BigInt::from(2));
    assert_eq!(matrix.b(), &BigInt::one());
    assert_eq!(matrix.c(), &BigInt::from(3));
    assert_eq!(matrix.d(), &BigInt::from(2));
}

#[test]
fn compose_with_identity_preserves_matrix() {
    let matrix = ModularMatrix::new(2, 1, 3, 2).expect("matrix should be valid");

    assert_eq!(
        matrix.compose(&ModularMatrix::identity()),
        Ok(matrix.clone())
    );
    assert_eq!(ModularMatrix::identity().compose(&matrix), Ok(matrix));
}

#[test]
fn s_squared_is_negative_identity() {
    let negative_identity = ModularMatrix::new(-1, 0, 0, -1).expect("matrix should be valid");

    assert_eq!(
        ModularMatrix::s().compose(&ModularMatrix::s()),
        Ok(negative_identity)
    );
}

#[test]
fn inverse_undoes_the_matrix_on_both_sides() {
    let matrix = ModularMatrix::new(2, 1, 3, 2).expect("matrix should be valid");
    let inverse = matrix.inverse().expect("inverse should stay in range");

    assert_eq!(matrix.compose(&inverse), Ok(ModularMatrix::identity()));
    assert_eq!(inverse.compose(&matrix), Ok(ModularMatrix::identity()));
}

#[test]
fn unbounded_entries_compose_and_invert_exactly() {
    let huge_shift = BigInt::one() << 200usize;
    let left = ModularMatrix::new(1, huge_shift.clone(), 0, 1).expect("matrix should be valid");
    let right = ModularMatrix::new(1, -huge_shift, 0, 1).expect("matrix should be valid");

    assert_eq!(left.compose(&right), Ok(ModularMatrix::identity()));
    assert_eq!(left.inverse(), Ok(right));
}

#[test]
fn apply_rejects_entries_too_large_for_f64() {
    let huge_shift = BigInt::one() << 20_000usize;
    let matrix = ModularMatrix::new(1, huge_shift, 0, 1).expect("matrix should be valid");

    assert_eq!(
        matrix.apply(&UpperHalfPlanePoint::tau_i()),
        Err(AnalyticCurveError::NumericalComparisonFailed)
    );
}

#[test]
fn t_translates_tau_by_one() {
    let tau = UpperHalfPlanePoint::tau_i();
    let image = ModularMatrix::t()
        .apply(&tau)
        .expect("T should preserve the upper half-plane");

    assert_eq!(image.tau(), &Complex64::new(1.0, 1.0));
}

#[test]
fn s_fixes_tau_i() {
    let tau = UpperHalfPlanePoint::tau_i();
    let image = ModularMatrix::s()
        .apply(&tau)
        .expect("S should preserve the upper half-plane");

    assert!((image.real_part() - 0.0).abs() < 1.0e-12);
    assert!((image.imaginary_part() - 1.0).abs() < 1.0e-12);
}

#[test]
fn composed_action_matches_sequential_action() {
    let tau = UpperHalfPlanePoint::tau_generic_example();
    let s = ModularMatrix::s();
    let t = ModularMatrix::t();
    let composed = s.compose(&t).expect("composition should stay in range");

    let sequential = s
        .apply(
            &t.apply(&tau)
                .expect("T should preserve the upper half-plane"),
        )
        .expect("S should preserve the upper half-plane");
    let direct = composed
        .apply(&tau)
        .expect("composed action should preserve the upper half-plane");

    assert!((sequential.real_part() - direct.real_part()).abs() < 1.0e-12);
    assert!((sequential.imaginary_part() - direct.imaginary_part()).abs() < 1.0e-12);
}

#[test]
fn modular_invariance_report_for_s_at_tau_i_is_exactly_stable() {
    let report = ModularMatrix::s()
        .verify_j_invariance_at(
            UpperHalfPlanePoint::tau_i(),
            LatticeSumTruncation::larger_for_comparison(),
            ApproxTolerance::strict(),
        )
        .expect("report should be computable");

    assert_eq!(report.original_tau(), &UpperHalfPlanePoint::tau_i());
    assert_eq!(report.transformed_tau(), &UpperHalfPlanePoint::tau_i());
    assert_eq!(report.matrix(), ModularMatrix::s());
    assert_eq!(
        report.truncation(),
        LatticeSumTruncation::larger_for_comparison()
    );
    assert_eq!(report.tolerance(), ApproxTolerance::strict());
    assert_eq!(report.original_j(), report.transformed_j());
    assert_eq!(report.difference(), &Complex64::new(0.0, 0.0));
    assert_eq!(report.absolute_difference(), 0.0);
    assert!(report.invariant_approximately());
}

#[test]
fn modular_invariance_report_tracks_transformed_tau_and_finite_residuals() {
    let tau = UpperHalfPlanePoint::tau_generic_example();
    let matrix = ModularMatrix::t();
    let report = matrix
        .verify_j_invariance_at(
            tau.clone(),
            LatticeSumTruncation::larger_for_comparison(),
            ApproxTolerance::loose(),
        )
        .expect("report should be computable");

    assert_eq!(report.original_tau(), &tau);
    assert_eq!(
        report.transformed_tau(),
        &UpperHalfPlanePoint::from_re_im(tau.real_part() + 1.0, tau.imaginary_part()).unwrap()
    );
    assert_eq!(report.matrix(), matrix);
    assert!(report.original_j().re.is_finite());
    assert!(report.original_j().im.is_finite());
    assert!(report.transformed_j().re.is_finite());
    assert!(report.transformed_j().im.is_finite());
    assert!(report.absolute_difference().is_finite());
}

#[test]
fn j_is_invariant_under_t_approximately() {
    let report = ModularMatrix::t()
        .verify_j_invariance_at(
            UpperHalfPlanePoint::tau_i(),
            LatticeSumTruncation::new(12).unwrap(),
            ApproxTolerance::new(1.0e-3, 1.0e-3),
        )
        .unwrap();

    assert!(report.invariant_approximately());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn pseudo_random_sl2_matrices_always_have_determinant_one(
        matrix in modular_matrix_strategy(),
    ) {
        prop_assert_eq!(matrix.determinant(), BigInt::one());
    }

    #[test]
    fn inverse_recovers_the_identity_for_generated_matrices(
        matrix in modular_matrix_strategy(),
    ) {
        let inverse = matrix.inverse().expect("bounded generated matrices should invert in range");

        prop_assert_eq!(matrix.compose(&inverse), Ok(ModularMatrix::identity()));
        prop_assert_eq!(inverse.compose(&matrix), Ok(ModularMatrix::identity()));
    }

    #[test]
    fn composition_is_associative_on_bounded_generated_words(
        left in modular_matrix_strategy(),
        middle in modular_matrix_strategy(),
        right in modular_matrix_strategy(),
    ) {
        let left_then_middle = left.compose(&middle).expect("bounded products should stay in range");
        let middle_then_right = middle.compose(&right).expect("bounded products should stay in range");

        prop_assert_eq!(
            left_then_middle.compose(&right),
            left.compose(&middle_then_right),
        );
    }

    #[test]
    fn modular_action_keeps_tau_in_the_upper_half_plane(
        matrix in modular_matrix_strategy(),
        tau in arb_upper_half_plane_point(AnalyticStrategyConfig::default()),
    ) {
        let image = matrix.apply(&tau).expect("SL_2(ℤ) should preserve the upper half-plane");

        prop_assert!(image.imaginary_part() > 0.0);
    }

    #[test]
    fn inverse_action_recovers_generic_tau(
        matrix in modular_matrix_strategy(),
        tau in arb_upper_half_plane_point(AnalyticStrategyConfig::default()),
    ) {
        let image = matrix.apply(&tau).expect("SL_2(ℤ) should preserve the upper half-plane");
        let recovered = matrix
            .inverse()
            .expect("bounded generated matrices should invert in range")
            .apply(&image)
            .expect("inverse action should preserve the upper half-plane");

        prop_assert!((recovered.real_part() - tau.real_part()).abs() < 1.0e-10);
        prop_assert!((recovered.imaginary_part() - tau.imaginary_part()).abs() < 1.0e-10);
    }

    #[test]
    fn modular_invariance_report_stays_finite_for_bounded_generated_inputs(
        matrix in modular_matrix_strategy(),
        tau in arb_upper_half_plane_point(AnalyticStrategyConfig::default()),
    ) {
        let report = matrix.verify_j_invariance_at(tau, LatticeSumTruncation::larger_for_comparison(), ApproxTolerance::loose(),
        ).expect("report should be computable for bounded generated inputs");

        prop_assert!(report.original_j().re.is_finite());
        prop_assert!(report.original_j().im.is_finite());
        prop_assert!(report.transformed_j().re.is_finite());
        prop_assert!(report.transformed_j().im.is_finite());
        prop_assert!(report.absolute_difference().is_finite());
    }

    #[test]
    fn pseudo_random_and_word_generated_matrices_compose_to_another_valid_matrix(
        (gamma, delta) in modular_matrix_pair_strategy(),
    ) {
        let composed = gamma.compose(&delta).expect("bounded generated products should stay in range");
        prop_assert_eq!(composed.determinant(), BigInt::one());
    }
}
