use core::num::NonZeroU32;
use num_bigint::{BigInt, BigUint};
use proptest::prelude::*;

use crate::elliptic_curves::{
    endomorphisms::{
        EndomorphismRingCandidateSet, EndomorphismRingReport,
        quadratic_orders::{
            DiscriminantSign, ImaginaryQuadraticOrder, ImaginaryQuadraticOrderError,
            QuadraticDiscriminant, QuadraticDiscriminantFactorizationError,
            QuadraticDiscriminantMod4, QuadraticOrderIndexError,
        },
    },
    frobenius::{FrobeniusCharacteristicPolynomial, FrobeniusTrace},
    traits::EnumerableCurveModel,
};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;
use crate::numerics::PositivePrimeError;
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_endomorphism_report_case,
};

fn f43_trace_three_characteristic_polynomial() -> FrobeniusCharacteristicPolynomial {
    let base_field =
        FiniteFieldDescriptor::new(43, NonZeroU32::new(1).expect("1 should be non-zero"))
            .expect("F_43 metadata should be internally consistent");
    let trace = FrobeniusTrace::from_order(base_field, 41).expect("t = 3 should be valid over F43");
    trace.characteristic_polynomial()
}

#[test]
fn frobenius_constructor_recovers_t_squared_minus_four_q() {
    let discriminant = QuadraticDiscriminant::from_frobenius_trace_and_field_order(3, 43);

    assert_eq!(discriminant.value(), &BigInt::from(-163));
}

#[test]
fn frobenius_polynomial_constructor_matches_the_existing_polynomial_discriminant() {
    let polynomial = f43_trace_three_characteristic_polynomial();
    let discriminant = QuadraticDiscriminant::from_frobenius_characteristic_polynomial(&polynomial);

    assert_eq!(
        discriminant.value(),
        &BigInt::from(polynomial.discriminant())
    );
    assert_eq!(discriminant.value(), &BigInt::from(-163));
}

#[test]
fn sign_helpers_distinguish_negative_zero_and_positive() {
    let negative = QuadraticDiscriminant::new(-163);
    let zero = QuadraticDiscriminant::new(0);
    let positive = QuadraticDiscriminant::new(5);

    assert_eq!(negative.sign(), DiscriminantSign::Negative);
    assert!(negative.is_negative());
    assert!(!negative.is_zero());
    assert!(!negative.is_positive());

    assert_eq!(zero.sign(), DiscriminantSign::Zero);
    assert!(!zero.is_negative());
    assert!(zero.is_zero());
    assert!(!zero.is_positive());

    assert_eq!(positive.sign(), DiscriminantSign::Positive);
    assert!(!positive.is_negative());
    assert!(!positive.is_zero());
    assert!(positive.is_positive());
}

#[test]
fn mod_four_classification_distinguishes_zero_one_and_other_residues() {
    let minus_163 = QuadraticDiscriminant::new(-163);
    let eight = QuadraticDiscriminant::new(8);
    let two = QuadraticDiscriminant::new(2);
    let minus_five = QuadraticDiscriminant::new(-5);

    assert_eq!(minus_163.mod_4_class(), QuadraticDiscriminantMod4::One);
    assert!(minus_163.is_congruent_to_1_mod_4());
    assert!(!minus_163.is_congruent_to_0_mod_4());

    assert_eq!(eight.mod_4_class(), QuadraticDiscriminantMod4::Zero);
    assert!(eight.is_congruent_to_0_mod_4());
    assert!(!eight.is_congruent_to_1_mod_4());

    assert_eq!(two.mod_4_class(), QuadraticDiscriminantMod4::Other(2));
    assert!(!two.is_congruent_to_0_mod_4());
    assert!(!two.is_congruent_to_1_mod_4());

    assert_eq!(
        minus_five.mod_4_class(),
        QuadraticDiscriminantMod4::Other(3)
    );
}

#[test]
fn fundamental_test_accepts_classical_negative_and_positive_examples() {
    assert!(QuadraticDiscriminant::new(-163).is_fundamental());
    assert!(QuadraticDiscriminant::new(-4).is_fundamental());
    assert!(QuadraticDiscriminant::new(5).is_fundamental());
    assert!(QuadraticDiscriminant::new(8).is_fundamental());
    assert!(QuadraticDiscriminant::new(12).is_fundamental());
}

#[test]
fn fundamental_test_rejects_non_fundamental_or_wrong_mod_four_inputs() {
    assert!(!QuadraticDiscriminant::new(0).is_fundamental());
    assert!(!QuadraticDiscriminant::new(1).is_fundamental());
    assert!(!QuadraticDiscriminant::new(2).is_fundamental());
    assert!(!QuadraticDiscriminant::new(3).is_fundamental());
    assert!(!QuadraticDiscriminant::new(16).is_fundamental());
    assert!(!QuadraticDiscriminant::new(-48).is_fundamental());
}

#[test]
fn big_field_orders_are_supported_without_i128_overflow_checks() {
    let discriminant = QuadraticDiscriminant::from_frobenius_trace_and_field_order(0, u128::MAX);

    assert!(discriminant.is_negative());
    assert!(discriminant.is_congruent_to_0_mod_4());
}

#[test]
fn factorization_of_a_fundamental_negative_discriminant_is_trivial() {
    let factorization = QuadraticDiscriminant::new(-163)
        .factorization()
        .expect("-163 should factor as a fundamental imaginary quadratic discriminant");

    assert_eq!(
        factorization.discriminant(),
        &QuadraticDiscriminant::new(-163)
    );
    assert_eq!(factorization.conductor(), &num_bigint::BigUint::from(1u8));
    assert_eq!(
        factorization.fundamental_discriminant(),
        &QuadraticDiscriminant::new(-163)
    );
    assert!(factorization.is_fundamental_already());
}

#[test]
fn factorization_extracts_the_square_part_in_a_non_fundamental_example() {
    let factorization = QuadraticDiscriminant::new(-48)
        .factorization()
        .expect("-48 should decompose as v^2 times a fundamental discriminant");

    assert_eq!(factorization.conductor(), &num_bigint::BigUint::from(4u8));
    assert_eq!(
        factorization.fundamental_discriminant(),
        &QuadraticDiscriminant::new(-3)
    );
    assert!(!factorization.is_fundamental_already());
}

#[test]
fn factorization_handles_the_even_branch_with_fundamental_discriminant_minus_four() {
    let factorization = QuadraticDiscriminant::new(-16)
        .factorization()
        .expect("-16 should decompose through the D_K = -4 branch");

    assert_eq!(factorization.conductor(), &num_bigint::BigUint::from(2u8));
    assert_eq!(
        factorization.fundamental_discriminant(),
        &QuadraticDiscriminant::new(-4)
    );
}

#[test]
fn factorization_can_build_the_maximal_order() {
    let factorization = QuadraticDiscriminant::new(-48)
        .factorization()
        .expect("-48 should admit a canonical factorization");
    let maximal_order = factorization
        .maximal_order()
        .expect("the factorization should recover the maximal order O_K");

    assert_eq!(
        maximal_order.fundamental_discriminant(),
        &QuadraticDiscriminant::new(-3)
    );
    assert_eq!(maximal_order.conductor(), &num_bigint::BigUint::from(1u8));
    assert_eq!(
        maximal_order.discriminant(),
        &QuadraticDiscriminant::new(-3)
    );
    assert!(maximal_order.is_maximal());
}

#[test]
fn factorization_rejects_zero_and_positive_inputs_for_the_imaginary_case() {
    assert_eq!(
        QuadraticDiscriminant::new(0).factorization(),
        Err(QuadraticDiscriminantFactorizationError::ZeroDiscriminant)
    );
    assert_eq!(
        QuadraticDiscriminant::new(12).factorization(),
        Err(QuadraticDiscriminantFactorizationError::PositiveDiscriminant)
    );
}

#[test]
fn factorization_rejects_negative_integers_that_are_not_quadratic_order_discriminants() {
    assert_eq!(
        QuadraticDiscriminant::new(-2).factorization(),
        Err(QuadraticDiscriminantFactorizationError::InvalidQuadraticOrderDiscriminant)
    );
}

#[test]
fn factorization_error_display_is_mathematically_specific() {
    assert_eq!(
        QuadraticDiscriminantFactorizationError::ZeroDiscriminant.to_string(),
        "quadratic discriminant factorization is undefined for Δ = 0"
    );
    assert_eq!(
        QuadraticDiscriminantFactorizationError::PositiveDiscriminant.to_string(),
        "quadratic discriminant factorization currently expects an imaginary discriminant Δ < 0"
    );
    assert_eq!(
        QuadraticDiscriminantFactorizationError::InvalidQuadraticOrderDiscriminant.to_string(),
        "quadratic discriminant is not congruent to 0 or 1 modulo 4, so it does not define a quadratic order"
    );
}

#[test]
fn imaginary_quadratic_order_new_builds_the_expected_discriminant() {
    let order = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(4u8),
    )
    .expect("negative fundamental discriminant with positive conductor should define an order");

    assert_eq!(
        order.fundamental_discriminant(),
        &QuadraticDiscriminant::new(-3)
    );
    assert_eq!(order.conductor(), &num_bigint::BigUint::from(4u8));
    assert_eq!(order.discriminant(), &QuadraticDiscriminant::new(-48));
    assert!(!order.is_maximal());
}

#[test]
fn imaginary_quadratic_order_from_factorization_roundtrips_the_existing_split() {
    let factorization = QuadraticDiscriminant::new(-48)
        .factorization()
        .expect("-48 should admit a canonical factorization");
    let order = ImaginaryQuadraticOrder::from_factorization(factorization)
        .expect("the canonical factorization of -48 should define an order");

    assert_eq!(
        order.fundamental_discriminant(),
        &QuadraticDiscriminant::new(-3)
    );
    assert_eq!(order.conductor(), &num_bigint::BigUint::from(4u8));
    assert_eq!(order.discriminant(), &QuadraticDiscriminant::new(-48));
}

#[test]
fn imaginary_quadratic_order_from_discriminant_uses_the_factorization_pipeline() {
    let order = ImaginaryQuadraticOrder::from_discriminant(&QuadraticDiscriminant::new(-16))
        .expect("-16 should define the order with conductor 2 inside D_K = -4");

    assert_eq!(
        order.fundamental_discriminant(),
        &QuadraticDiscriminant::new(-4)
    );
    assert_eq!(order.conductor(), &num_bigint::BigUint::from(2u8));
    assert_eq!(order.discriminant(), &QuadraticDiscriminant::new(-16));
}

#[test]
fn imaginary_quadratic_order_new_rejects_invalid_inputs() {
    assert_eq!(
        ImaginaryQuadraticOrder::new(
            QuadraticDiscriminant::new(5),
            num_bigint::BigUint::from(1u8)
        ),
        Err(ImaginaryQuadraticOrderError::NonNegativeFundamentalDiscriminant)
    );
    assert_eq!(
        ImaginaryQuadraticOrder::new(
            QuadraticDiscriminant::new(-48),
            num_bigint::BigUint::from(1u8)
        ),
        Err(ImaginaryQuadraticOrderError::NonFundamentalDiscriminant)
    );
    assert_eq!(
        ImaginaryQuadraticOrder::new(
            QuadraticDiscriminant::new(-3),
            num_bigint::BigUint::from(0u8)
        ),
        Err(ImaginaryQuadraticOrderError::ZeroConductor)
    );
}

#[test]
fn imaginary_quadratic_order_from_discriminant_rejects_non_imaginary_cases() {
    assert_eq!(
        ImaginaryQuadraticOrder::from_discriminant(&QuadraticDiscriminant::new(12)),
        Err(ImaginaryQuadraticOrderError::NonImaginaryOrderDiscriminant)
    );
}

#[test]
fn imaginary_quadratic_order_error_display_is_mathematically_specific() {
    assert_eq!(
        ImaginaryQuadraticOrderError::NonNegativeFundamentalDiscriminant.to_string(),
        "imaginary quadratic orders require a negative fundamental discriminant D_K < 0"
    );
    assert_eq!(
        ImaginaryQuadraticOrderError::NonFundamentalDiscriminant.to_string(),
        "imaginary quadratic order construction requires a fundamental discriminant D_K"
    );
    assert_eq!(
        ImaginaryQuadraticOrderError::ZeroConductor.to_string(),
        "imaginary quadratic order construction requires a positive conductor f >= 1"
    );
    assert_eq!(
        ImaginaryQuadraticOrderError::NonImaginaryOrderDiscriminant.to_string(),
        "quadratic discriminant does not define an imaginary quadratic order"
    );
    assert_eq!(
        QuadraticOrderIndexError::DifferentQuadraticFields.to_string(),
        "quadratic-order index is defined only for orders in the same imaginary quadratic field"
    );
    assert_eq!(
        QuadraticOrderIndexError::NotSuborder.to_string(),
        "quadratic-order index requires an inclusion O_f ⊆ O_g"
    );
}

#[test]
fn imaginary_quadratic_orders_can_compare_field_membership_and_containment() {
    let maximal = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(1u8),
    )
    .expect("maximal order should construct");
    let conductor_two = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(2u8),
    )
    .expect("conductor-two order should construct");
    let conductor_four = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(4u8),
    )
    .expect("conductor-four order should construct");
    let other_field = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-4),
        num_bigint::BigUint::from(1u8),
    )
    .expect("other maximal order should construct");

    assert!(maximal.same_quadratic_field(&conductor_two));
    assert!(!maximal.same_quadratic_field(&other_field));

    assert!(conductor_four.is_suborder_of(&conductor_two));
    assert!(conductor_four.is_suborder_of(&maximal));
    assert!(maximal.is_overorder_of(&conductor_four));
    assert!(conductor_two.is_overorder_of(&conductor_four));

    assert!(!conductor_two.is_suborder_of(&conductor_four));
    assert!(!conductor_two.is_suborder_of(&other_field));
}

#[test]
fn imaginary_quadratic_order_index_matches_the_conductor_quotient() {
    let maximal = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(1u8),
    )
    .expect("maximal order should construct");
    let conductor_two = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(2u8),
    )
    .expect("conductor-two order should construct");
    let conductor_four = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(4u8),
    )
    .expect("conductor-four order should construct");

    assert_eq!(
        maximal.index_of_suborder(&conductor_four),
        Ok(num_bigint::BigUint::from(4u8))
    );
    assert_eq!(
        conductor_two.index_of_suborder(&conductor_four),
        Ok(num_bigint::BigUint::from(2u8))
    );
}

#[test]
fn imaginary_quadratic_order_index_rejects_non_inclusions_and_other_fields() {
    let maximal = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(1u8),
    )
    .expect("maximal order should construct");
    let conductor_two = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-3),
        num_bigint::BigUint::from(2u8),
    )
    .expect("conductor-two order should construct");
    let other_field = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-4),
        num_bigint::BigUint::from(1u8),
    )
    .expect("other maximal order should construct");

    assert_eq!(
        conductor_two.index_of_suborder(&maximal),
        Err(QuadraticOrderIndexError::NotSuborder)
    );
    assert_eq!(
        maximal.index_of_suborder(&other_field),
        Err(QuadraticOrderIndexError::DifferentQuadraticFields)
    );
}

#[test]
fn frobenius_generated_order_is_contained_in_the_maximal_order_from_the_same_factorization() {
    let factorization = QuadraticDiscriminant::new(-48)
        .factorization()
        .expect("-48 should admit a canonical factorization");
    let frobenius_order = ImaginaryQuadraticOrder::from_factorization(factorization.clone())
        .expect("the factorization should define Z[pi]");
    let maximal_order = factorization
        .maximal_order()
        .expect("the same factorization should recover O_K");

    assert!(frobenius_order.is_suborder_of(&maximal_order));
    assert!(maximal_order.is_overorder_of(&frobenius_order));
}

#[test]
fn candidate_set_enumerates_all_orders_between_frobenius_and_maximal() {
    let factorization = QuadraticDiscriminant::new(-48)
        .factorization()
        .expect("-48 should admit a canonical factorization");
    let candidates = factorization
        .endomorphism_ring_candidates()
        .expect("the factorization should enumerate every intermediate order");

    assert_eq!(
        candidates
            .candidate_orders()
            .iter()
            .map(|order| order.conductor().clone())
            .collect::<Vec<_>>(),
        vec![
            num_bigint::BigUint::from(1u8),
            num_bigint::BigUint::from(2u8),
            num_bigint::BigUint::from(4u8),
        ]
    );
    assert_eq!(candidates.len(), 3);
    assert!(candidates.maximal_order().is_maximal());
    assert_eq!(
        candidates.frobenius_order().discriminant(),
        &QuadraticDiscriminant::new(-48)
    );
    assert!(
        candidates
            .frobenius_order()
            .is_suborder_of(candidates.maximal_order())
    );
}

#[test]
fn candidate_set_from_discriminant_uses_the_same_factorization_data() {
    let candidates =
        EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(-16))
            .expect("-16 should produce the conductor-divisor candidate orders");

    assert_eq!(
        candidates.frobenius_conductor(),
        &num_bigint::BigUint::from(2u8)
    );
    assert_eq!(
        candidates.fundamental_discriminant(),
        &QuadraticDiscriminant::new(-4)
    );
    assert_eq!(candidates.candidate_orders().len(), 2);
    assert_eq!(
        candidates.maximal_order().discriminant(),
        &QuadraticDiscriminant::new(-4)
    );
    assert_eq!(
        candidates.frobenius_order().discriminant(),
        &QuadraticDiscriminant::new(-16)
    );
}

#[test]
fn candidate_set_can_measure_indices_over_the_frobenius_order() {
    let candidates =
        EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(-48))
            .expect("-48 should produce the conductor-divisor candidate orders");

    assert_eq!(
        candidates.index_over_frobenius_order(candidates.maximal_order()),
        Ok(num_bigint::BigUint::from(4u8))
    );
    assert_eq!(
        candidates.index_over_frobenius_order(&candidates.candidate_orders()[1]),
        Ok(num_bigint::BigUint::from(2u8))
    );
    assert_eq!(
        candidates.index_over_frobenius_order(candidates.frobenius_order()),
        Ok(num_bigint::BigUint::from(1u8))
    );
}

#[test]
fn factorization_local_view_records_the_prime_adic_conductor_gap() {
    let factorization = QuadraticDiscriminant::new(-144)
        .factorization()
        .expect("-144 should factor as 6^2 times -4");

    let local_two = factorization
        .local_view_at(&BigUint::from(2u8))
        .expect("2 should define a valid local view");
    let local_three = factorization
        .local_view_at(&BigUint::from(3u8))
        .expect("3 should define a valid local view");
    let local_five = factorization
        .local_view_at(&BigUint::from(5u8))
        .expect("5 should define a valid local view");

    assert_eq!(local_two.prime(), &BigUint::from(2u8));
    assert_eq!(local_two.frobenius_conductor_valuation(), 1);
    assert!(!local_two.is_trivial());
    assert_eq!(
        local_two.local_conductor_exponents().collect::<Vec<_>>(),
        vec![0, 1]
    );

    assert_eq!(local_three.frobenius_conductor_valuation(), 1);
    assert!(!local_three.is_trivial());

    assert_eq!(local_five.frobenius_conductor_valuation(), 0);
    assert!(local_five.is_trivial());
    assert_eq!(
        local_five.local_conductor_exponents().collect::<Vec<_>>(),
        vec![0]
    );
}

#[test]
fn candidate_set_local_view_delegates_to_the_same_factorization_data() {
    let candidates =
        EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(-144))
            .expect("-144 should produce Frobenius-compatible candidate orders");

    let local_view = candidates
        .local_view_at(&BigUint::from(2u8))
        .expect("2 should define a valid local view");

    assert_eq!(local_view.frobenius_conductor_valuation(), 1);
    assert_eq!(
        local_view.local_conductor_exponents().collect::<Vec<_>>(),
        vec![0, 1]
    );
}

#[test]
fn local_view_rejects_non_prime_inputs() {
    let factorization = QuadraticDiscriminant::new(-48)
        .factorization()
        .expect("-48 should admit a canonical factorization");

    assert_eq!(
        factorization.local_view_at(&BigUint::from(0u8)),
        Err(PositivePrimeError::Zero)
    );
    assert_eq!(
        factorization.local_view_at(&BigUint::from(1u8)),
        Err(PositivePrimeError::One)
    );
    assert_eq!(
        factorization.local_view_at(&BigUint::from(4u8)),
        Err(PositivePrimeError::Composite)
    );
}

#[test]
fn imaginary_quadratic_order_volcanic_level_matches_the_local_conductor_valuation() {
    let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-3), BigUint::from(12u8))
        .expect("conductor-twelve order should construct");

    let level_two = order
        .volcanic_level_at(&BigUint::from(2u8))
        .expect("2 should define a valid volcanic level candidate");
    let level_three = order
        .volcanic_level_at(&BigUint::from(3u8))
        .expect("3 should define a valid volcanic level candidate");
    let level_five = order
        .volcanic_level_at(&BigUint::from(5u8))
        .expect("5 should define a valid volcanic level candidate");

    assert_eq!(level_two.prime(), &BigUint::from(2u8));
    assert_eq!(level_two.order(), &order);
    assert_eq!(level_two.level(), 2);
    assert_eq!(level_two.conductor(), &BigUint::from(12u8));
    assert_eq!(level_two.discriminant(), &QuadraticDiscriminant::new(-432));

    assert_eq!(level_three.level(), 1);
    assert_eq!(level_five.level(), 0);
}

#[test]
fn candidate_set_can_annotate_each_order_with_its_volcanic_level_candidate() {
    let candidates =
        EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(-144))
            .expect("-144 should produce Frobenius-compatible candidate orders");

    let level_candidates = candidates
        .volcanic_level_candidates_at(&BigUint::from(2u8))
        .expect("2 should define valid volcanic level candidates");

    let levels: Vec<_> = level_candidates
        .iter()
        .map(|candidate| (candidate.conductor().clone(), candidate.level()))
        .collect();

    assert_eq!(
        levels,
        vec![
            (BigUint::from(1u8), 0),
            (BigUint::from(2u8), 1),
            (BigUint::from(3u8), 0),
            (BigUint::from(6u8), 1),
        ]
    );
}

#[test]
fn volcanic_level_candidates_reject_non_prime_inputs() {
    let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-3), BigUint::from(4u8))
        .expect("conductor-four order should construct");

    assert_eq!(
        order.volcanic_level_at(&BigUint::from(0u8)),
        Err(PositivePrimeError::Zero)
    );
    assert_eq!(
        order.volcanic_level_at(&BigUint::from(1u8)),
        Err(PositivePrimeError::One)
    );
    assert_eq!(
        order.volcanic_level_at(&BigUint::from(4u8)),
        Err(PositivePrimeError::Composite)
    );
}

#[test]
fn candidate_set_exposes_hasse_cover_relations_for_a_branching_example() {
    let candidates =
        EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(-144))
            .expect("-144 should produce the conductor-divisor candidate orders");

    let edges: Vec<_> = candidates
        .hasse_cover_relations()
        .into_iter()
        .map(|edge| {
            (
                edge.overorder().conductor().clone(),
                edge.suborder().conductor().clone(),
                edge.index().clone(),
            )
        })
        .collect();

    assert_eq!(
        edges,
        vec![
            (
                num_bigint::BigUint::from(1u8),
                num_bigint::BigUint::from(2u8),
                num_bigint::BigUint::from(2u8),
            ),
            (
                num_bigint::BigUint::from(1u8),
                num_bigint::BigUint::from(3u8),
                num_bigint::BigUint::from(3u8),
            ),
            (
                num_bigint::BigUint::from(2u8),
                num_bigint::BigUint::from(6u8),
                num_bigint::BigUint::from(3u8),
            ),
            (
                num_bigint::BigUint::from(3u8),
                num_bigint::BigUint::from(6u8),
                num_bigint::BigUint::from(2u8),
            ),
        ]
    );
}

#[test]
fn candidate_set_index_over_frobenius_order_rejects_orders_from_other_fields() {
    let candidates =
        EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(-48))
            .expect("-48 should produce the conductor-divisor candidate orders");
    let other_field = ImaginaryQuadraticOrder::new(
        QuadraticDiscriminant::new(-4),
        num_bigint::BigUint::from(1u8),
    )
    .expect("other maximal order should construct");

    assert_eq!(
        candidates.index_over_frobenius_order(&other_field),
        Err(QuadraticOrderIndexError::DifferentQuadraticFields)
    );
}

#[test]
fn endomorphism_ring_report_packages_only_frobenius_compatible_candidates() {
    let base_field =
        FiniteFieldDescriptor::new(43, NonZeroU32::new(1).expect("1 should be non-zero"))
            .expect("F_43 metadata should be internally consistent");
    let trace = FrobeniusTrace::from_order(base_field, 41).expect("t = 3 should be valid over F43");
    let report = EndomorphismRingReport::from_frobenius_trace(trace)
        .expect("ordinary Frobenius data should define a compatible candidate report");

    assert!(matches!(
        report,
        EndomorphismRingReport::OrdinaryQuadraticOrderCandidates { .. }
    ));
    assert_eq!(
        report
            .factorization()
            .expect("ordinary branch should expose a factorization")
            .discriminant(),
        &QuadraticDiscriminant::new(-163)
    );
    assert_eq!(report.candidate_count(), Some(1));
    assert_eq!(
        report.candidate_orders(),
        report.candidate_set().map(|set| set.candidate_orders())
    );
    assert_eq!(report.frobenius_order(), report.maximal_order());
    assert_eq!(report.sandwich_inclusion_holds(), Some(true));
    assert!(report.is_ordinary());
    assert!(!report.is_supersingular());
}

#[test]
fn endomorphism_ring_report_handles_nontrivial_intermediate_candidates() {
    let base_field =
        FiniteFieldDescriptor::new(13, NonZeroU32::new(1).expect("1 should be non-zero"))
            .expect("F_13 metadata should be internally consistent");
    let trace =
        FrobeniusTrace::from_order(base_field, 18).expect("t = -4 should be valid over F13");
    let report = EndomorphismRingReport::from_frobenius_trace(trace)
        .expect("negative discriminant should produce a compatible candidate report");

    assert!(matches!(
        report,
        EndomorphismRingReport::OrdinaryQuadraticOrderCandidates { .. }
    ));
    assert_eq!(
        report
            .candidate_orders()
            .expect("ordinary branch should expose candidate orders")
            .iter()
            .map(|order| order.conductor().clone())
            .collect::<Vec<_>>(),
        vec![
            num_bigint::BigUint::from(1u8),
            num_bigint::BigUint::from(3u8),
        ]
    );
    assert_eq!(
        report
            .frobenius_order()
            .expect("ordinary branch should expose Z[pi]")
            .discriminant(),
        &QuadraticDiscriminant::new(-36)
    );
    assert_eq!(
        report
            .maximal_order()
            .expect("ordinary branch should expose O_K")
            .discriminant(),
        &QuadraticDiscriminant::new(-4)
    );
    assert_eq!(report.sandwich_inclusion_holds(), Some(true));
}

#[test]
fn endomorphism_ring_report_distinguishes_the_supersingular_branch() {
    let base_field =
        FiniteFieldDescriptor::new(5, NonZeroU32::new(1).expect("1 should be non-zero"))
            .expect("F_5 metadata should be internally consistent");
    let trace = FrobeniusTrace::from_order(base_field, 6).expect("t = 0 should be valid over F5");
    let report = EndomorphismRingReport::from_frobenius_trace(trace)
        .expect("supersingular Frobenius data should still produce a report");

    assert!(matches!(
        report,
        EndomorphismRingReport::SupersingularQuaternionicPlaceholder { .. }
    ));
    assert_eq!(report.frobenius_discriminant().trace(), BigInt::from(0));
    assert_eq!(report.factorization(), None);
    assert_eq!(report.candidate_set(), None);
    assert!(report.is_supersingular());
    assert!(!report.is_ordinary());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn generated_endomorphism_reports_stay_frobenius_compatible(
        case in arb_endomorphism_report_case::<crate::fields::Fp17>(CurveStrategyConfig::default()),
    ) {
        let discriminant = case.report.frobenius_discriminant();
        let trace = discriminant.frobenius_trace();

        prop_assert_eq!(trace.curve_order(), BigUint::from(case.curve.order() as u64));
        prop_assert_eq!(discriminant.trace(), trace.trace());
        prop_assert_eq!(case.report.is_ordinary(), !case.report.is_supersingular());

        if case.report.is_ordinary() {
            prop_assert!(case.report.factorization().is_some());
            prop_assert!(case.report.candidate_set().is_some());
            prop_assert_eq!(case.report.sandwich_inclusion_holds(), Some(true));
        } else {
            prop_assert!(case.report.factorization().is_none());
            prop_assert!(case.report.candidate_set().is_none());
        }
    }
}
