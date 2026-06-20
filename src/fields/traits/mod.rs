//! Trait surfaces for field-like structures and field capabilities.
//!
//! The public trait namespace intentionally lives under
//! `crate::fields::traits::...` so the top-level `fields` module remains a
//! barrel of concrete mathematical families rather than a mixed bag of both
//! families and abstractions.

pub mod ambient_field;
pub mod artin_schreier;
pub mod cbrt_field;
pub mod enumerative_finite_field;
pub mod field;
pub mod finite_field;
pub mod pth_root_extraction;
pub mod quadratic_character;
pub mod sqrt_field;

pub use ambient_field::AmbientField;
pub use artin_schreier::CharacteristicTwoArtinSchreierField;
pub use cbrt_field::CbrtField;
pub use enumerative_finite_field::EnumerableFiniteField;
pub use field::Field;
pub use finite_field::FiniteField;
pub use pth_root_extraction::PthRootExtraction;
pub use quadratic_character::{QuadraticCharacterFiniteField, QuadraticCharacterValue};
pub use sqrt_field::SqrtField;
