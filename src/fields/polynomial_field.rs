use core::marker::PhantomData;

use crate::DensePolynomial;
use crate::fields::{errors::FieldError, traits::Field};
use crate::polynomials::{
    IrreducibilityBackend, IrreducibilityStatus, PolynomialError, irreducibility_status,
};

type DenseTriple<F> = (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>);

/// Modulus polynomial used to define a quotient of `F[x]`.
pub struct PolynomialModulus<F: Field> {
    coefficients: Vec<F::Elem>,
    field: PhantomData<F>,
}

impl<F: Field> Clone for PolynomialModulus<F> {
    fn clone(&self) -> Self {
        Self {
            coefficients: self.coefficients.clone(),
            field: PhantomData,
        }
    }
}

impl<F: Field> core::fmt::Debug for PolynomialModulus<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PolynomialModulus")
            .field("coefficients", &self.coefficients)
            .finish()
    }
}

impl<F: Field> PolynomialModulus<F> {
    /// Creates a modulus polynomial from coefficients in ascending order.
    ///
    /// This constructor currently performs only the weakest structural check:
    /// the modulus must be non-constant.
    ///
    /// That is enough to talk about a quotient algebra `F[x] / (m(x))`, but it
    /// is not enough to guarantee that the quotient behaves as a field.
    ///
    /// For quotient-field constructions, the modulus polynomial is typically
    /// expected to be irreducible over the coefficient field `F`. That stronger
    /// requirement is intentionally not enforced by the constructor itself,
    /// because:
    ///
    /// - this type is also useful for generic quotient-algebra scaffolding
    /// - some base-field backends only support partial irreducibility checks
    /// - the project still wants to keep the low-level constructor lightweight
    ///
    /// In other words:
    ///
    /// - `new(...)` means “this is a syntactically valid non-constant modulus”
    /// - `check_field_modulus_requirements(...)` is the stronger question
    ///   “is this suitable for a quotient field in the currently supported
    ///   backend?”
    pub fn new(coefficients: Vec<F::Elem>) -> Result<Self, FieldError> {
        if coefficients.len() < 2 {
            return Err(FieldError::InvalidPolynomialModulus);
        }

        Ok(Self {
            coefficients,
            field: PhantomData,
        })
    }

    /// Returns the raw coefficients.
    pub fn coefficients(&self) -> &[F::Elem] {
        &self.coefficients
    }

    /// Returns the degree of the modulus polynomial.
    pub fn degree(&self) -> usize {
        self.coefficients.len().saturating_sub(1)
    }
}

impl<F: Field + IrreducibilityBackend> PolynomialModulus<F> {
    /// Checks whether the modulus satisfies the stronger conditions usually
    /// required for quotient-field constructions.
    ///
    /// The current implementation delegates to the dense polynomial
    /// irreducibility infrastructure in `polynomials` and translates its
    /// structured result back into the field-domain error surface.
    ///
    /// Current behavior:
    ///
    /// - `Constant` should be unreachable because [`PolynomialModulus::new`]
    ///   already rejects constant inputs
    /// - `Linear` and `Irreducible` are accepted
    /// - any reducibility witness, or a theoretical reducibility conclusion
    ///   without witness, is rejected as
    ///   [`FieldError::NonIrreduciblePolynomial`]
    /// - exact-but-partial backends such as `Q` may report
    ///   [`FieldError::UndeterminedPolynomialModulusIrreducibility`]
    pub fn check_field_modulus_requirements(&self) -> Result<(), FieldError> {
        let dense_modulus = DensePolynomial::<F>::new(self.coefficients.clone());

        match irreducibility_status(&dense_modulus) {
            Ok(IrreducibilityStatus::Linear | IrreducibilityStatus::Irreducible) => Ok(()),
            Ok(
                IrreducibilityStatus::Reducible { .. }
                | IrreducibilityStatus::ReducibleWithoutWitness { .. },
            ) => Err(FieldError::NonIrreduciblePolynomial),
            Ok(IrreducibilityStatus::Constant) => Err(FieldError::InvalidPolynomialModulus),
            Err(PolynomialError::UndeterminedIrreducibility(_)) => {
                Err(FieldError::UndeterminedPolynomialModulusIrreducibility)
            }
            Err(PolynomialError::InvalidBaseField(field_error)) => Err(field_error),
            Err(PolynomialError::UnsupportedIrreducibilityBackend(message)) => {
                Err(FieldError::Unsupported(message))
            }
            Err(_) => Err(FieldError::Unsupported(
                "unexpected polynomial-domain error during modulus irreducibility checking",
            )),
        }
    }
}

/// Element of a quotient algebra `F[x] / (m(x))`.
///
/// Mathematically, this type is meant to model a polynomial representative
/// together with the modulus polynomial that defines the quotient relation.
/// Two representatives should be understood as denoting the same quotient
/// element when they differ by a multiple of `m(x)`.
///
/// In other words, if `m(x)` is the chosen modulus polynomial, then
/// polynomials `f(x)` and `g(x)` represent the same quotient element whenever:
///
/// `f(x) - g(x) = q(x) * m(x)`
///
/// for some polynomial `q(x)` over the same coefficient field `F`.
///
/// This struct deliberately stores:
///
/// - a raw representative polynomial through `coefficients`
/// - the modulus polynomial that defines the quotient space
///
/// The current representation uses coefficients in ascending degree order.
/// For example:
///
/// - `[a0, a1, a2]` represents `a0 + a1*x + a2*x^2`
/// - a modulus `[b0, b1, b2]` represents `b0 + b1*x + b2*x^2`
///
/// This type is now an operational quotient-element layer, but it remains more
/// lightweight and pedagogical than [`crate::fields::ExtensionField`]:
///
/// - representatives are allowed to be stored unreduced and can later be
///   canonicalized through [`PolynomialFieldElement::reduce`]
/// - equality is interpreted by quotient class when the defining modulus
///   matches
/// - addition, subtraction, negation, multiplication, inversion, and division
///   are available directly on the value type
/// - irreducibility is still checked explicitly, not enforced by the
///   constructor
///
/// That means `PolynomialFieldElement<F>` should currently be read as
/// “an autocontained quotient value with explicit representative and opt-in
/// field validation” rather than as the primary field-family backend.
pub struct PolynomialFieldElement<F: Field> {
    coefficients: Vec<F::Elem>,
    modulus: PolynomialModulus<F>,
}

impl<F: Field> PolynomialFieldElement<F> {
    /// Builds a quotient element from a polynomial representative.
    ///
    /// The `coefficients` slice stores the chosen representative in ascending
    /// degree order, while `modulus` defines the quotient relation.
    ///
    /// At the moment this constructor validates only the basic structural
    /// requirement that the modulus polynomial is non-constant. It does not yet:
    ///
    /// - reduce the representative modulo the modulus
    /// - remove trailing zero coefficients
    /// - verify that the modulus is irreducible
    ///
    /// Those behaviors are intentionally split across the current scaffold:
    ///
    /// - irreducibility can now be checked explicitly through
    ///   [`PolynomialFieldElement::check_field_conditions`] when the base field
    ///   implements the polynomial irreducibility backend
    /// - representative reduction is available through
    ///   [`PolynomialFieldElement::reduce`]
    pub fn new(
        coefficients: Vec<F::Elem>,
        modulus: PolynomialModulus<F>,
    ) -> Result<Self, FieldError> {
        if modulus.coefficients().is_empty() {
            return Err(FieldError::InvalidPolynomialModulus);
        }

        Ok(Self {
            coefficients,
            modulus,
        })
    }

    /// Returns the representative coefficients.
    ///
    /// The coefficients are stored in ascending degree order. The returned
    /// slice exposes the current explicit representative, which may still be
    /// unreduced until [`PolynomialFieldElement::reduce`] is called.
    pub fn coefficients(&self) -> &[F::Elem] {
        &self.coefficients
    }

    /// Returns the degree of the currently stored representative, if it is
    /// non-zero.
    pub fn degree(&self) -> Option<usize> {
        self.coefficients.len().checked_sub(1)
    }

    /// Returns whether the currently stored representative is the zero
    /// polynomial.
    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }

    /// Returns the defining modulus.
    ///
    /// Conceptually, this is the polynomial `m(x)` in the quotient
    /// construction `F[x] / (m(x))`.
    pub fn modulus(&self) -> &PolynomialModulus<F> {
        &self.modulus
    }

    /// Reduces the representative polynomial modulo the field relation.
    ///
    /// The current implementation computes the Euclidean remainder of the
    /// stored representative polynomial modulo the defining polynomial.
    ///
    /// This gives a canonical representative of degree strictly smaller than
    /// the modulus degree whenever the modulus is non-zero. It does not, by
    /// itself, certify that the quotient is a true field; that stronger
    /// question still belongs to
    /// [`PolynomialFieldElement::check_field_conditions`].
    pub fn reduce(&self) -> Result<Self, FieldError> {
        let representative = DensePolynomial::<F>::new(self.coefficients.clone());
        let modulus = DensePolynomial::<F>::new(self.modulus.coefficients().to_vec());
        let remainder = representative
            .rem(&modulus)
            .map_err(Self::map_polynomial_error)?;

        Self::new(remainder.coefficients().to_vec(), self.modulus.clone())
    }

    /// Returns a canonical remainder representative.
    ///
    /// This is a convenience alias for [`PolynomialFieldElement::reduce`] to
    /// make call sites read more naturally when the goal is to obtain a
    /// reduced quotient value.
    pub fn reduced(&self) -> Result<Self, FieldError> {
        self.reduce()
    }

    /// Returns whether the currently stored representative is already reduced
    /// modulo the defining polynomial.
    pub fn is_reduced(&self) -> Result<bool, FieldError> {
        let reduced = self.reduce()?;
        Ok(Self::same_coefficients(
            self.coefficients(),
            reduced.coefficients(),
        ))
    }

    /// Adds two quotient representatives and reduces the result canonically.
    pub fn add(&self, rhs: &Self) -> Result<Self, FieldError> {
        self.ensure_compatible_modulus(rhs)?;

        let sum = self.dense_representative().add(&rhs.dense_representative());

        Self::new(sum.coefficients().to_vec(), self.modulus.clone())?.reduce()
    }

    /// Negates the stored representative coefficient-wise and reduces the
    /// result canonically.
    pub fn neg(&self) -> Result<Self, FieldError> {
        let negated = DensePolynomial::<F>::new(self.coefficients.iter().map(F::neg).collect());
        Self::new(negated.coefficients().to_vec(), self.modulus.clone())?.reduce()
    }

    /// Subtracts two quotient representatives and reduces the result
    /// canonically.
    pub fn sub(&self, rhs: &Self) -> Result<Self, FieldError> {
        self.ensure_compatible_modulus(rhs)?;
        self.add(&rhs.neg()?)
    }

    /// Multiplies two quotient representatives and reduces the result modulo
    /// the defining polynomial.
    pub fn mul(&self, rhs: &Self) -> Result<Self, FieldError> {
        self.ensure_compatible_modulus(rhs)?;

        let product = self.dense_representative().mul(&rhs.dense_representative());

        Self::new(product.coefficients().to_vec(), self.modulus.clone())?.reduce()
    }

    /// Computes the multiplicative inverse when the element is invertible in
    /// the quotient.
    ///
    /// In a general quotient algebra, not every non-zero representative need
    /// be invertible. When the modulus is irreducible, every non-zero element
    /// becomes invertible and this method behaves like field inversion.
    pub fn inverse(&self) -> Result<Self, FieldError> {
        let reduced = self.reduce()?;
        if reduced.is_zero() {
            return Err(FieldError::DivisionByZero);
        }

        let modulus = self.dense_modulus();
        let representative = reduced.dense_representative();
        let (gcd, bezout, _) =
            Self::extended_gcd(representative, modulus).map_err(Self::map_polynomial_error)?;

        let Some(unit) = gcd.constant_term() else {
            return Err(FieldError::NonInvertibleElement);
        };

        if gcd.degree().is_some_and(|degree| degree > 0) {
            return Err(FieldError::NonInvertibleElement);
        }

        let unit_inverse = F::inverse(unit)?;
        let inverse = bezout.scale(&unit_inverse);

        Self::new(inverse.coefficients().to_vec(), self.modulus.clone())?.reduce()
    }

    /// Divides by another quotient representative when the divisor is
    /// invertible.
    pub fn div(&self, rhs: &Self) -> Result<Self, FieldError> {
        self.ensure_compatible_modulus(rhs)?;
        self.mul(&rhs.inverse()?)
    }

    /// Checks the modulus for the extra conditions required by a field.
    ///
    /// In a true quotient field, the modulus polynomial must satisfy stronger
    /// conditions than in a generic quotient algebra, most notably
    /// irreducibility over the coefficient field.
    ///
    /// This method is available only when the base field implements the
    /// polynomial irreducibility backend used by the crate.
    pub fn check_field_conditions(&self) -> Result<(), FieldError>
    where
        F: IrreducibilityBackend,
    {
        self.modulus.check_field_modulus_requirements()
    }

    fn ensure_compatible_modulus(&self, rhs: &Self) -> Result<(), FieldError> {
        if self.modulus == rhs.modulus {
            Ok(())
        } else {
            Err(FieldError::IncompatibleFieldParameters)
        }
    }

    fn dense_representative(&self) -> DensePolynomial<F> {
        DensePolynomial::<F>::new(self.coefficients.clone())
    }

    fn dense_modulus(&self) -> DensePolynomial<F> {
        DensePolynomial::<F>::new(self.modulus.coefficients().to_vec())
    }

    fn same_coefficients(lhs: &[F::Elem], rhs: &[F::Elem]) -> bool {
        lhs.len() == rhs.len() && lhs.iter().zip(rhs).all(|(a, b)| F::eq(a, b))
    }

    fn extended_gcd(
        left: DensePolynomial<F>,
        right: DensePolynomial<F>,
    ) -> Result<DenseTriple<F>, PolynomialError> {
        let mut old_r = left;
        let mut r = right;
        let mut old_s = DensePolynomial::<F>::constant(F::one());
        let mut s = DensePolynomial::<F>::new(Vec::new());
        let mut old_t = DensePolynomial::<F>::new(Vec::new());
        let mut t = DensePolynomial::<F>::constant(F::one());

        while !r.is_zero() {
            let (quotient, remainder) = old_r.div_rem(&r)?;
            old_r = r;
            r = remainder;

            let next_s = old_s.sub(&quotient.mul(&s));
            old_s = s;
            s = next_s;

            let next_t = old_t.sub(&quotient.mul(&t));
            old_t = t;
            t = next_t;
        }

        Ok((old_r, old_s, old_t))
    }

    fn map_polynomial_error(error: PolynomialError) -> FieldError {
        match error {
            PolynomialError::DivisionByZeroPolynomial => FieldError::InvalidPolynomialModulus,
            PolynomialError::NonInvertibleLeadingCoefficient => FieldError::NonInvertibleElement,
            _ => FieldError::Unsupported(
                "unexpected polynomial-domain error during quotient reduction",
            ),
        }
    }
}

impl<F: Field> Clone for PolynomialFieldElement<F> {
    fn clone(&self) -> Self {
        Self {
            coefficients: self.coefficients.clone(),
            modulus: self.modulus.clone(),
        }
    }
}

impl<F: Field> core::fmt::Debug for PolynomialFieldElement<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PolynomialFieldElement")
            .field("coefficients", &self.coefficients)
            .field("modulus", &self.modulus)
            .finish()
    }
}

impl<F: Field> PartialEq for PolynomialModulus<F> {
    fn eq(&self, other: &Self) -> bool {
        PolynomialFieldElement::<F>::same_coefficients(self.coefficients(), other.coefficients())
    }
}

impl<F: Field> PartialEq for PolynomialFieldElement<F> {
    fn eq(&self, other: &Self) -> bool {
        if self.modulus != other.modulus {
            return false;
        }

        match (self.reduce(), other.reduce()) {
            (Ok(lhs), Ok(rhs)) => Self::same_coefficients(lhs.coefficients(), rhs.coefficients()),
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::{PolynomialFieldElement, PolynomialModulus};
    use crate::fields::{ComplexApprox, Field, FieldError, Fp, Q};
    use num_bigint::BigInt;
    use num_rational::BigRational;

    type F17 = Fp<17>;

    #[test]
    fn quotient_modulus_requires_degree_at_least_one() {
        let error = PolynomialModulus::<F17>::new(vec![F17::elem_from_u64(42)])
            .expect_err("constant modulus should fail");
        assert!(matches!(
            error,
            crate::fields::FieldError::InvalidPolynomialModulus
        ));
    }

    #[test]
    fn quotient_modulus_preserves_coefficients_and_degree() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(5),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");

        let coefficients = modulus.coefficients();
        assert_eq!(coefficients.len(), 4);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(3)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(0)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(5)));
        assert!(F17::eq(&coefficients[3], &F17::elem_from_u64(1)));
        assert_eq!(modulus.degree(), 3);
    }

    #[test]
    fn quotient_element_preserves_representative_coefficients() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");

        let element = PolynomialFieldElement::<F17>::new(
            vec![
                F17::elem_from_u64(9),
                F17::elem_from_u64(4),
                F17::elem_from_u64(15),
            ],
            modulus,
        )
        .expect("element should exist");

        let coefficients = element.coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(9)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(4)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(15)));
    }

    #[test]
    fn quotient_element_exposes_defining_modulus() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(2),
            F17::elem_from_u64(3),
            F17::elem_from_u64(4),
        ])
        .expect("modulus should exist");

        let element = PolynomialFieldElement::<F17>::new(vec![F17::elem_from_u64(8)], modulus)
            .expect("element should exist");

        let stored_modulus = element.modulus();
        assert_eq!(stored_modulus.degree(), 2);
        let coefficients = stored_modulus.coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(2)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(3)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(4)));
    }

    #[test]
    fn quotient_element_can_store_zero_representative() {
        let modulus =
            PolynomialModulus::<F17>::new(vec![F17::elem_from_u64(1), F17::elem_from_u64(1)])
                .expect("modulus should exist");

        let element = PolynomialFieldElement::<F17>::new(Vec::new(), modulus)
            .expect("zero representative should be allowed");

        assert!(element.coefficients().is_empty());
    }

    #[test]
    fn quotient_modulus_field_requirements_accept_irreducible_prime_field_modulus() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        modulus
            .check_field_modulus_requirements()
            .expect("x^2 + 3 is irreducible over F17");
    }

    #[test]
    fn quotient_modulus_field_requirements_reject_reducible_prime_field_modulus() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");

        assert_eq!(
            modulus.check_field_modulus_requirements(),
            Err(FieldError::NonIrreduciblePolynomial)
        );
    }

    #[test]
    fn quotient_element_field_conditions_use_real_irreducibility_checks() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let element = PolynomialFieldElement::<F17>::new(vec![F17::elem_from_u64(5)], modulus)
            .expect("element should exist");
        element
            .check_field_conditions()
            .expect("element should inherit irreducible modulus");
    }

    #[test]
    fn complex_modulus_requirements_reject_degree_two_polynomials() {
        let modulus = PolynomialModulus::<ComplexApprox>::new(vec![
            ComplexApprox::from_i64(1),
            ComplexApprox::zero(),
            ComplexApprox::one(),
        ])
        .expect("modulus should exist");

        assert_eq!(
            modulus.check_field_modulus_requirements(),
            Err(FieldError::NonIrreduciblePolynomial)
        );
    }

    #[test]
    fn rational_modulus_requirements_can_succeed_on_certified_examples() {
        let modulus =
            PolynomialModulus::<Q>::new(vec![Q::from_i64(1), Q::from_i64(0), Q::from_i64(1)])
                .expect("modulus should exist");

        modulus
            .check_field_modulus_requirements()
            .expect("x^2 + 1 is irreducible over Q");
    }

    #[test]
    fn rational_modulus_requirements_report_inconclusive_cases_honestly() {
        let leading = [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
            .into_iter()
            .fold(BigInt::from(1_u64), |accumulator, prime| {
                accumulator * BigInt::from(prime)
            });
        let constant = &leading + BigInt::from(1_u64);
        let modulus = PolynomialModulus::<Q>::new(vec![
            BigRational::from_integer(constant),
            Q::zero(),
            Q::zero(),
            Q::zero(),
            BigRational::from_integer(leading),
        ])
        .expect("modulus should exist");

        assert_eq!(
            modulus.check_field_modulus_requirements(),
            Err(FieldError::UndeterminedPolynomialModulusIrreducibility)
        );
    }

    #[test]
    fn quotient_reduction_computes_canonical_remainder_over_prime_fields() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let element = PolynomialFieldElement::<F17>::new(
            vec![
                F17::elem_from_u64(1),
                F17::elem_from_u64(2),
                F17::elem_from_u64(3),
            ],
            modulus,
        )
        .expect("element should exist");
        let reduced = element.reduce().expect("reduction should succeed");

        let coefficients = reduced.coefficients();
        assert_eq!(coefficients.len(), 2);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(15)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(2)));
    }

    #[test]
    fn quotient_reduction_normalizes_to_zero_when_representative_is_a_multiple_of_modulus() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let element = PolynomialFieldElement::<F17>::new(
            vec![
                F17::elem_from_u64(2),
                F17::elem_from_u64(0),
                F17::elem_from_u64(2),
            ],
            modulus,
        )
        .expect("element should exist");

        let reduced = element.reduce().expect("reduction should succeed");
        assert!(reduced.coefficients().is_empty());
    }

    #[test]
    fn quotient_reduction_preserves_already_reduced_representatives() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let element = PolynomialFieldElement::<F17>::new(
            vec![F17::elem_from_u64(4), F17::elem_from_u64(7)],
            modulus,
        )
        .expect("element should exist");

        let reduced = element.reduce().expect("reduction should succeed");
        let coefficients = reduced.coefficients();
        assert_eq!(coefficients.len(), 2);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(4)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(7)));
    }

    #[test]
    fn quotient_element_reports_reduced_status_and_degree() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let unreduced = PolynomialFieldElement::<F17>::new(
            vec![
                F17::elem_from_u64(1),
                F17::elem_from_u64(2),
                F17::elem_from_u64(3),
            ],
            modulus.clone(),
        )
        .expect("element should exist");
        let reduced = unreduced.reduced().expect("reduction should succeed");

        assert_eq!(unreduced.degree(), Some(2));
        assert!(!unreduced.is_zero());
        assert!(!unreduced.is_reduced().expect("check should succeed"));

        assert_eq!(reduced.degree(), Some(1));
        assert!(reduced.is_reduced().expect("check should succeed"));
    }

    #[test]
    fn quotient_arithmetic_add_sub_neg_and_mul_reduce_canonically() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let left = PolynomialFieldElement::<F17>::new(
            vec![F17::elem_from_u64(1), F17::elem_from_u64(1)],
            modulus.clone(),
        )
        .expect("left should exist");
        let right = PolynomialFieldElement::<F17>::new(
            vec![F17::elem_from_u64(3), F17::elem_from_u64(16)],
            modulus,
        )
        .expect("right should exist");

        let sum = left.add(&right).expect("addition should succeed");
        let diff = left.sub(&right).expect("subtraction should succeed");
        let neg = right.neg().expect("negation should succeed");
        let product = left.mul(&right).expect("multiplication should succeed");

        assert!(PolynomialFieldElement::<F17>::same_coefficients(
            sum.coefficients(),
            &[F17::elem_from_u64(4)]
        ));
        assert!(PolynomialFieldElement::<F17>::same_coefficients(
            diff.coefficients(),
            &[F17::elem_from_u64(15), F17::elem_from_u64(2)]
        ));
        assert!(PolynomialFieldElement::<F17>::same_coefficients(
            neg.coefficients(),
            &[F17::elem_from_u64(14), F17::elem_from_u64(1)]
        ));
        assert!(PolynomialFieldElement::<F17>::same_coefficients(
            product.coefficients(),
            &[F17::elem_from_u64(4), F17::elem_from_u64(2)]
        ));
    }

    #[test]
    fn quotient_equality_is_by_reduced_class_not_raw_storage() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let x_squared = PolynomialFieldElement::<F17>::new(
            vec![
                F17::elem_from_u64(0),
                F17::elem_from_u64(0),
                F17::elem_from_u64(1),
            ],
            modulus.clone(),
        )
        .expect("element should exist");
        let minus_one = PolynomialFieldElement::<F17>::new(vec![F17::elem_from_u64(16)], modulus)
            .expect("element should exist");

        assert!(x_squared == minus_one);
    }

    #[test]
    fn quotient_operations_reject_incompatible_moduli() {
        let lhs = PolynomialFieldElement::<F17>::new(
            vec![F17::elem_from_u64(1)],
            PolynomialModulus::<F17>::new(vec![
                F17::elem_from_u64(1),
                F17::elem_from_u64(0),
                F17::elem_from_u64(1),
            ])
            .expect("lhs modulus should exist"),
        )
        .expect("lhs should exist");
        let rhs = PolynomialFieldElement::<F17>::new(
            vec![F17::elem_from_u64(1)],
            PolynomialModulus::<F17>::new(vec![
                F17::elem_from_u64(3),
                F17::elem_from_u64(0),
                F17::elem_from_u64(1),
            ])
            .expect("rhs modulus should exist"),
        )
        .expect("rhs should exist");

        assert_eq!(lhs.add(&rhs), Err(FieldError::IncompatibleFieldParameters));
        assert_eq!(lhs.mul(&rhs), Err(FieldError::IncompatibleFieldParameters));
        assert!(!(lhs == rhs));
    }

    #[test]
    fn quotient_inverse_and_division_work_when_modulus_is_irreducible() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let element = PolynomialFieldElement::<F17>::new(
            vec![F17::elem_from_u64(1), F17::elem_from_u64(1)],
            modulus.clone(),
        )
        .expect("element should exist");
        let inverse = element.inverse().expect("element should be invertible");
        let one = element.mul(&inverse).expect("product should succeed");
        let quotient = element.div(&element).expect("division should succeed");

        assert!(PolynomialFieldElement::<F17>::same_coefficients(
            one.coefficients(),
            &[F17::elem_from_u64(1)]
        ));
        assert!(PolynomialFieldElement::<F17>::same_coefficients(
            quotient.coefficients(),
            &[F17::elem_from_u64(1)]
        ));
    }

    #[test]
    fn quotient_inverse_rejects_non_units_in_reducible_quotients() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let zero_divisor = PolynomialFieldElement::<F17>::new(
            vec![F17::elem_from_u64(13), F17::elem_from_u64(1)],
            modulus,
        )
        .expect("element should exist");

        assert_eq!(
            zero_divisor.inverse(),
            Err(FieldError::NonInvertibleElement)
        );
    }

    fn irreducible_f17_modulus() -> PolynomialModulus<F17> {
        PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("x^2 + 3 should be a valid modulus")
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(32))]

        #[test]
        fn property_polynomial_field_reduction_is_idempotent(
            coefficients in prop::collection::vec(0u64..17, 0..=5),
        ) {
            let element = PolynomialFieldElement::<F17>::new(
                coefficients.into_iter().map(F17::elem_from_u64).collect(),
                irreducible_f17_modulus(),
            )
            .expect("element should exist");

            let reduced_once = element.reduce().expect("reduction should succeed");
            let reduced_twice = reduced_once.reduce().expect("reduction should succeed");

            prop_assert_eq!(reduced_once.clone(), reduced_twice);
            prop_assert!(reduced_once.is_reduced().expect("reduced check should succeed"));
        }

        #[test]
        fn property_polynomial_field_nonzero_elements_are_invertible(
            coefficients in prop::collection::vec(0u64..17, 0..=2),
        ) {
            let element = PolynomialFieldElement::<F17>::new(
                coefficients.into_iter().map(F17::elem_from_u64).collect(),
                irreducible_f17_modulus(),
            )
            .expect("element should exist");
            let reduced = element.reduce().expect("reduction should succeed");

            prop_assume!(!reduced.is_zero());

            let inverse = reduced.inverse().expect("non-zero element should be invertible over an irreducible modulus");
            let one = reduced.mul(&inverse).expect("product should succeed");

            prop_assert!(PolynomialFieldElement::<F17>::same_coefficients(
                one.coefficients(),
                &[F17::one()]
            ));
        }
    }
}
