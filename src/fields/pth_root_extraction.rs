use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::fields::{
    ExtensionField, ExtensionFieldElement, ExtensionFieldSpec, FiniteField, Fp, FpElem,
};

/// Capability trait for algebraic objects that can admit a `p`-th root in
/// characteristic `p`.
///
/// This trait is intentionally value-oriented rather than field-family-
/// oriented. That keeps the same surface usable for several mathematically
/// different kinds of objects, for example:
///
/// - elements of finite fields
/// - polynomials over finite fields
///
/// Contract:
///
/// - `pth_root(x)` returns `Some(y)` only when `y^p = x` in the relevant
///   ambient algebraic structure
/// - if no such `y` exists in that structure, it returns `None`
/// - when the `p`-th root is unique, implementations should return that root
pub trait PthRootExtraction: Sized {
    /// Returns one `p`-th root of `self` when it exists.
    fn pth_root(&self) -> Option<Self>;

    /// Returns whether `self` admits a `p`-th root in the same ambient
    /// algebraic structure.
    fn has_pth_root(&self) -> bool {
        self.pth_root().is_some()
    }
}

/// Returns the exact exponent `p^(n-1)` for `F_(p^n)`.
pub fn finite_field_pth_root_exponent<F: FiniteField>() -> BigUint {
    BigUint::from(F::characteristic()).pow(F::extension_degree().get().saturating_sub(1))
}

/// Raises a finite-field element to an arbitrary-precision unsigned exponent
/// by repeated squaring.
pub fn finite_field_pow_biguint<F: FiniteField>(x: &F::Elem, exponent: &BigUint) -> F::Elem {
    let mut result = F::one();
    let mut base = x.clone();
    let mut exp = exponent.clone();

    while !exp.is_zero() {
        if (&exp & BigUint::one()) == BigUint::one() {
            result = F::mul(&result, &base);
        }

        exp >>= 1usize;

        if !exp.is_zero() {
            base = F::square(&base);
        }
    }

    result
}

/// Applies the characteristic-Frobenius map `x ↦ x^p` in a finite field.
pub fn finite_field_frobenius_p<F: FiniteField>(x: &F::Elem) -> F::Elem {
    F::pow(x, F::characteristic())
}

/// Returns the unique `p`-th root of a finite-field element.
///
/// If `F = F_(p^n)`, then `Fr_p^{-1}(x) = x^(p^(n-1))`.
pub fn finite_field_pth_root<F: FiniteField>(x: &F::Elem) -> F::Elem {
    finite_field_pow_biguint::<F>(x, &finite_field_pth_root_exponent::<F>())
}

impl<const P: u64> PthRootExtraction for FpElem<P> {
    fn pth_root(&self) -> Option<Self> {
        Some(finite_field_pth_root::<Fp<P>>(self))
    }
}

impl<S: ExtensionFieldSpec> PthRootExtraction for ExtensionFieldElement<S>
where
    ExtensionField<S>: FiniteField<Elem = ExtensionFieldElement<S>>,
{
    fn pth_root(&self) -> Option<Self> {
        Some(finite_field_pth_root::<ExtensionField<S>>(self))
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use crate::fields::{
        EnumerableFiniteField, Field, Fp, PthRootExtraction, finite_field_frobenius_p,
        finite_field_pow_biguint, finite_field_pth_root_exponent,
    };

    type F17 = Fp<17>;

    crate::fields::define_fp_quadratic_extension!(
        spec: F17Sqrt3PthRootSpec,
        field: F17Sqrt3PthRoot,
        base: F17,
        non_residue: 3,
        name: "F17(sqrt(3)) for p-th-root tests",
    );

    #[test]
    fn finite_field_prime_elements_have_identity_pth_root() {
        for element in F17::elements() {
            assert_eq!(element.pth_root(), Some(element));
            assert!(element.has_pth_root());
        }
    }

    #[test]
    fn finite_field_frobenius_matches_the_pth_power_map() {
        for element in F17Sqrt3PthRoot::elements() {
            assert_eq!(
                finite_field_frobenius_p::<F17Sqrt3PthRoot>(&element),
                F17Sqrt3PthRoot::pow(&element, F17Sqrt3PthRoot::characteristic())
            );
        }
    }

    #[test]
    fn finite_field_elements_admit_unique_pth_roots() {
        for element in F17Sqrt3PthRoot::elements() {
            let root = element
                .pth_root()
                .expect("finite-field elements should always admit a p-th root");

            assert_eq!(finite_field_frobenius_p::<F17Sqrt3PthRoot>(&root), element);
            assert!(element.has_pth_root());
        }
    }

    #[test]
    fn finite_field_pth_root_uses_the_biguint_exponent_formula() {
        let exponent = finite_field_pth_root_exponent::<F17Sqrt3PthRoot>();

        assert_eq!(exponent, BigUint::from(17u8));

        for element in F17Sqrt3PthRoot::elements() {
            assert_eq!(
                element.pth_root(),
                Some(finite_field_pow_biguint::<F17Sqrt3PthRoot>(
                    &element, &exponent
                ))
            );
        }
    }
}
