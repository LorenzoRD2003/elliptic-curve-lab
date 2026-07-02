/// Defines an educational quadratic extension over a prime-field family.
///
/// This helper is meant for small examples, tests, and walkthrough code that
/// wants a concrete type-level presentation of
/// `Fp<M, LIMBS>[x] / (x^2 - d)`
/// without rewriting the same [`ExtensionFieldSpec`] boilerplate.
///
/// The generated specification validates the modulus through
/// [`PolynomialModulus::check_field_modulus_requirements`], so the caller is
/// still responsible for choosing a value of `d` that is genuinely
/// non-square when a true quadratic field extension is desired.
///
/// Example:
///
/// ```ignore
/// use elliptic_algorithms_lab::fields::{traits::Field};
///
/// type F19 = crate::fields::Fp19;
///
/// elliptic_algorithms_lab::fields::extension_field::define_fp_quadratic_extension!(
///     spec: F19Sqrt2Spec,
///     field: F19Sqrt2,
///     base: F19,
///     non_residue: 2,
///     name: "F19(sqrt(2))",
/// );
/// ```
#[doc(hidden)]
#[macro_export]
macro_rules! define_fp_quadratic_extension {
    (
        spec: $spec:ident,
        field: $field:ident,
        base: $base:ty,
        non_residue: $non_residue:expr,
        name: $name:expr $(,)?
    ) => {
        #[derive(Clone, Copy)]
        struct $spec;

        impl $crate::fields::extension_field::ExtensionFieldSpec for $spec {
            type Base = $base;

            const NAME: &'static str = $name;

            fn defining_modulus() -> $crate::fields::polynomial_field::PolynomialModulus<Self::Base>
            {
                $crate::fields::polynomial_field::PolynomialModulus::<Self::Base>::new(vec![
                    <Self::Base as $crate::fields::traits::Field>::from_i64(-($non_residue)),
                    <Self::Base as $crate::fields::traits::Field>::zero(),
                    <Self::Base as $crate::fields::traits::Field>::one(),
                ])
                .expect("x^2 - d should be a valid structural modulus")
            }

            fn check_field_conditions() -> Result<(), $crate::fields::FieldError> {
                Self::defining_modulus().check_field_modulus_requirements()
            }
        }

        type $field = $crate::fields::extension_field::ExtensionField<$spec>;
    };
}

/// Defines an educational quadratic extension over `Q`.
///
/// This helper is meant for small exact examples such as `Q(sqrt(2))` that
/// would otherwise repeat the same [`ExtensionFieldSpec`] boilerplate.
///
/// The generated specification validates the modulus through
/// [`PolynomialModulus::check_field_modulus_requirements`], so it remains
/// honest about whether `x^2 - d` really defines a field extension.
///
/// Example:
///
/// ```ignore
/// elliptic_algorithms_lab::fields::extension_field::define_q_quadratic_extension!(
///     spec: QSqrt2Spec,
///     field: QSqrt2,
///     radicand: 2,
///     name: "Q(sqrt(2))",
/// );
/// ```
#[doc(hidden)]
#[macro_export]
macro_rules! define_q_quadratic_extension {
    (
        spec: $spec:ident,
        field: $field:ident,
        radicand: $radicand:expr,
        name: $name:expr $(,)?
    ) => {
        #[derive(Clone, Copy)]
        struct $spec;

        impl $crate::fields::extension_field::ExtensionFieldSpec for $spec {
            type Base = $crate::fields::Q;

            const NAME: &'static str = $name;

            fn defining_modulus() -> $crate::fields::polynomial_field::PolynomialModulus<Self::Base>
            {
                $crate::fields::polynomial_field::PolynomialModulus::<Self::Base>::new(vec![
                    <$crate::fields::Q as $crate::fields::traits::Field>::from_i64(-($radicand)),
                    <$crate::fields::Q as $crate::fields::traits::Field>::zero(),
                    <$crate::fields::Q as $crate::fields::traits::Field>::one(),
                ])
                .expect("x^2 - d should be a valid structural modulus")
            }

            fn check_field_conditions() -> Result<(), $crate::fields::FieldError> {
                Self::defining_modulus().check_field_modulus_requirements()
            }
        }

        type $field = $crate::fields::extension_field::ExtensionField<$spec>;
    };
}
