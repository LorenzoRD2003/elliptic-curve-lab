//! Foundational scaffolding for mathematical and cryptographic algorithms.
//!
//! The crate intentionally starts with small, documented interfaces and
//! lightweight placeholder implementations so the core abstractions can evolve
//! with tests before the heavy algebraic algorithms arrive.

// pub mod algorithms;
pub mod elliptic_curves;
pub mod fields;
pub mod polynomials;
pub mod visualization;

pub use elliptic_curves::{
    AffineCurveModel, AffinePoint, CurveError, CurveModel, EnumerableCurveModel,
    FiniteAbelianGroupStructure, FiniteGroupCurveModel, GroupCurveModel, LiftXCoordinate,
    PointIndexSampler, ShortWeierstrassCurve,
};
pub use fields::{
    ComplexApprox, EnumerableFiniteField, ExtensionField, ExtensionFieldElement,
    ExtensionFieldSpec, Field, FieldError, FiniteField, FiniteFieldDescriptor, Fp, FpElem,
    PolynomialFieldElement, PolynomialModulus, Q, SqrtField, addition_table, describe_complex,
    describe_complex_polynomial_modulus_as_field_modulus, describe_extension_field,
    describe_extension_field_element, describe_prime_polynomial_field_element,
    describe_prime_polynomial_modulus, describe_prime_polynomial_modulus_as_field_modulus,
    describe_rational, explain_add, explain_complex_polynomial_modulus_irreducibility,
    explain_complex_square_root, explain_extension_field_add, explain_extension_field_inverse,
    explain_extension_field_mul, explain_extension_field_reduction, explain_inverse, explain_mul,
    explain_prime_field_square_root, explain_prime_polynomial_field_add,
    explain_prime_polynomial_field_inverse, explain_prime_polynomial_field_mul,
    explain_prime_polynomial_field_reduction, explain_prime_polynomial_modulus_irreducibility,
    explain_prime_polynomial_storage, explain_rational_add, explain_rational_div,
    explain_rational_inverse, explain_rational_mul, explain_rational_square_root, format_complex,
    format_complex_polynomial, format_extension_field, format_extension_field_element,
    format_fp_elem, format_prime_field, format_prime_polynomial,
    format_prime_polynomial_field_element, format_prime_polynomial_modulus, format_rational,
    format_rational_field, inverses_table, multiplication_table,
};
pub use polynomials::{
    DensePolynomial, IrreducibilityBackend, IrreducibilityStatus, PolynomialError,
    ReducibilityReason, SparsePolynomial, VisualizablePolynomial, describe_dense_polynomial,
    describe_multivariate_polynomial, describe_sparse_polynomial, explain_dense_storage,
    explain_multivariate_storage, explain_sparse_storage, format_dense_polynomial, format_monomial,
    format_multivariate_polynomial, format_sparse_polynomial, irreducibility_status,
    is_irreducible,
};
pub use visualization::{Visualizable, VisualizableField};
pub use visualization::{
    describe_curve, describe_group_structure, describe_membership, describe_order_distribution,
    describe_point, describe_point_order, describe_scalar_mul, explain_add as explain_curve_add,
    explain_point_order, format_curve, format_point, format_point_compact, list_points,
    summarize_group_structure, summarize_order_distribution,
};
