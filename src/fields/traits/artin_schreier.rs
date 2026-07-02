use crate::fields::{
    error::FieldError,
    traits::{EnumerableFiniteField, Field, FiniteField},
};

type ArtinSchreierSolutionPair<F> = (<F as Field>::Elem, <F as Field>::Elem);
type ArtinSchreierPairResult<F> = Result<Option<ArtinSchreierSolutionPair<F>>, FieldError>;

/// Capability trait for solving characteristic-`2` Artin-Schreier equations
///
/// `z^2 + z = a`
///
/// over small finite fields that can honestly enumerate all their elements.
///
/// This trait is intentionally separate from [`super::SqrtField`]:
///
/// - square roots solve equations of the form `z^2 = a`
/// - Artin-Schreier equations solve the affine-linearized form `z^2 + z = a`
///
/// In characteristic `2`, the latter is the right normalization for many
/// general-Weierstrass `x`-fibers when the coefficient of `y` is non-zero.
///
/// The current default algorithms are intentionally educational:
///
/// - trace computation uses the defining sum `a + a^2 + ... + a^{2^{n-1}}`
/// - existence uses that trace criterion
/// - witness recovery uses exhaustive search after the criterion says a
///   solution should exist
pub trait CharacteristicTwoArtinSchreierField: FiniteField + EnumerableFiniteField {
    /// Returns the absolute trace
    ///
    /// `Tr_{F_{2^n}/F_2}(a) = a + a^2 + ... + a^{2^{n-1}}`
    ///
    /// as an element of the ambient field.
    ///
    /// For characteristic-`2` finite fields this value lands in the prime
    /// subfield `{0, 1}` embedded into the represented field family.
    ///
    /// Complexity: `Θ(n)` field additions and squarings, where
    /// `n = [F : F_2]` is the extension degree of the finite field.
    fn artin_schreier_trace_to_prime_field(a: &Self::Elem) -> Result<Self::Elem, FieldError> {
        Self::ensure_characteristic_two()?;

        let mut trace = Self::zero();
        let mut term = a.clone();
        for _ in 0..Self::extension_degree().get() {
            trace = Self::add(&trace, &term);
            term = Self::square(&term);
        }

        Ok(trace)
    }

    /// Returns whether `z^2 + z = a` has a solution in the field.
    ///
    /// Over `F_{2^n}`, this happens exactly when the absolute trace of `a` is
    /// zero.
    ///
    /// Complexity: `Θ(n)` field additions and squarings under the current
    /// trace-based criterion.
    fn artin_schreier_has_solution(a: &Self::Elem) -> Result<bool, FieldError> {
        Ok(Self::is_zero(&Self::artin_schreier_trace_to_prime_field(
            a,
        )?))
    }

    /// Returns one solution to `z^2 + z = a` when it exists.
    ///
    /// The current educational implementation first uses the trace criterion
    /// to rule out the impossible case, then finds a witness by exhaustive
    /// search over the finite field.
    ///
    /// Complexity:
    /// - unsolvable case: `Θ(n)` field additions and squarings
    /// - solvable case: `Θ(n + q)` field additions, squarings, and equality
    ///   checks in the worst case, where `q = #F`
    fn solve_artin_schreier(a: &Self::Elem) -> Result<Option<Self::Elem>, FieldError> {
        Self::ensure_characteristic_two()?;

        if !Self::artin_schreier_has_solution(a)? {
            return Ok(None);
        }

        let solution = Self::elements()
            .into_iter()
            .find(|candidate| Self::eq(&Self::add(&Self::square(candidate), candidate), a));

        solution.map(Some).ok_or(FieldError::Unsupported(
            "Artin-Schreier trace criterion predicted a solution, but exhaustive search did not find one",
        ))
    }

    /// Returns the two solutions to `z^2 + z = a` when they exist.
    ///
    /// In characteristic `2`, if `z` is one solution then `z + 1` is the
    /// other one.
    ///
    /// Complexity: the cost of [`Self::solve_artin_schreier`] plus `Θ(1)` one
    /// extra field addition to produce the second solution.
    fn solve_artin_schreier_pair(a: &Self::Elem) -> ArtinSchreierPairResult<Self> {
        let solution = match Self::solve_artin_schreier(a)? {
            Some(solution) => solution,
            None => return Ok(None),
        };

        let partner = Self::add(&solution, &Self::one());
        Ok(Some((solution, partner)))
    }

    fn ensure_characteristic_two() -> Result<(), FieldError> {
        Self::check_structure()?;
        if !Self::has_characteristic(2) {
            return Err(FieldError::Unsupported(
                "Artin-Schreier solving is only implemented for finite fields of characteristic 2",
            ));
        }
        Ok(())
    }
}

impl<F> CharacteristicTwoArtinSchreierField for F where F: FiniteField + EnumerableFiniteField {}

#[cfg(test)]
mod tests {
    use crate::fields::traits::*;

    use crate::fields::{
        extension_field::{ExtensionField, ExtensionFieldSpec},
        polynomial_field::PolynomialModulus,
        traits::CharacteristicTwoArtinSchreierField,
    };

    type F2 = crate::fields::Fp2;

    #[derive(Clone, Copy)]
    struct F4ArtinSchreierSpec;

    impl ExtensionFieldSpec for F4ArtinSchreierSpec {
        type Base = F2;

        const NAME: &'static str = "F4 for Artin-Schreier tests";

        fn defining_modulus() -> PolynomialModulus<Self::Base> {
            PolynomialModulus::<Self::Base>::new(vec![F2::one(), F2::one(), F2::one()])
                .expect("x^2 + x + 1 should be a valid structural modulus")
        }

        fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
            Self::defining_modulus().check_field_modulus_requirements()
        }
    }

    type F4 = ExtensionField<F4ArtinSchreierSpec>;

    #[test]
    fn prime_field_trace_is_the_identity_over_f2() {
        assert!(F2::eq(
            &F2::artin_schreier_trace_to_prime_field(&F2::zero()).expect("F2 should be supported"),
            &F2::zero()
        ));
        assert!(F2::eq(
            &F2::artin_schreier_trace_to_prime_field(&F2::one()).expect("F2 should be supported"),
            &F2::one()
        ));
    }

    #[test]
    fn prime_field_solver_finds_both_roots_of_z_squared_plus_z_equals_zero() {
        let (left, right) = F2::solve_artin_schreier_pair(&F2::zero())
            .expect("F2 should be supported")
            .expect("z^2 + z = 0 should have solutions in F2");

        assert!(!F2::eq(&left, &right));
        assert!(F2::eq(&F2::add(&F2::square(&left), &left), &F2::zero()));
        assert!(F2::eq(&F2::add(&F2::square(&right), &right), &F2::zero()));
    }

    #[test]
    fn prime_field_solver_rejects_z_squared_plus_z_equals_one_over_f2() {
        assert!(!F2::artin_schreier_has_solution(&F2::one()).expect("F2 should be supported"));
        assert_eq!(
            F2::solve_artin_schreier(&F2::one()).expect("F2 should be supported"),
            None
        );
    }

    #[test]
    fn extension_field_trace_detects_when_artin_schreier_has_a_solution() {
        let alpha = F4::element(vec![F2::zero(), F2::one()]);

        assert!(F4::eq(
            &F4::artin_schreier_trace_to_prime_field(&F4::one()).expect("F4 should be supported"),
            &F4::zero()
        ));
        assert!(F4::eq(
            &F4::artin_schreier_trace_to_prime_field(&alpha).expect("F4 should be supported"),
            &F4::one()
        ));

        assert!(F4::artin_schreier_has_solution(&F4::one()).expect("F4 should be supported"));
        assert!(!F4::artin_schreier_has_solution(&alpha).expect("F4 should be supported"));
    }

    #[test]
    fn extension_field_solver_returns_a_solution_pair_and_their_difference_is_one() {
        let (left, right) = F4::solve_artin_schreier_pair(&F4::one())
            .expect("F4 should be supported")
            .expect("z^2 + z = 1 should have solutions in F4");

        assert!(F4::eq(&F4::add(&F4::square(&left), &left), &F4::one()));
        assert!(F4::eq(&F4::add(&F4::square(&right), &right), &F4::one()));
        assert!(F4::eq(&F4::sub(&right, &left), &F4::one()));
    }

    #[test]
    fn non_characteristic_two_fields_are_rejected_honestly() {
        type F5 = crate::fields::Fp5;

        assert_eq!(
            F5::solve_artin_schreier(&F5::one()),
            Err(crate::fields::FieldError::Unsupported(
                "Artin-Schreier solving is only implemented for finite fields of characteristic 2",
            ))
        );
    }
}
