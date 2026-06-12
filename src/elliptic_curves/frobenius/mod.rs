//! The current module keeps two distinctions explicit:
//!
//! - the absolute Frobenius `π_p`, which raises coordinates to `p`-power
//! - the relative Frobenius `π_q`, where `q = p^r` is the size of the
//!   chosen finite base field
//!
//! For a curve defined over `F_q`, `π_q` is an endomorphism of that curve.
//! By contrast, `π_p` need only land on the `p`-power Frobenius twist in
//! general. When the curve coefficients already lie in the prime field, those
//! two curve models coincide, but the fixed points can still differ: this is
//! the first visible distinction between `E(F_p)` and larger finite-field
//! point sets represented in the same coordinate field.
//!
//! In the current educational implementation, both maps act only on point
//! coordinates already represented in one concrete finite field backend. The
//! module does not yet introduce a separate point type for geometric points
//! over an algebraic closure. Even so, the fixed-point story is already
//! visible:
//!
//! - points fixed by `π_p` behave like `E(F_p)` inside the chosen backend
//! - points fixed only by `π_q` can witness larger rationality fields such
//!   as `E(F_{p^r})`

mod character_sum;
mod characteristic_equation;
mod characteristic_polynomial;
mod curve_type;
mod discriminant;
mod extension_counts;
mod group_order;
mod hasse;
mod hasse_multiple_search;
mod interval;
mod isogeny;
mod orbit;
mod quadratic_twist;
mod short_weierstrass;
mod torsion;
mod torsion_matrix;
mod trace;
mod types;
mod zeta;

pub use crate::elliptic_curves::traits::FrobeniusTraceCurveModel;
pub use character_sum::CharacterSumPointCount;
pub use characteristic_equation::{
    FrobeniusCharacteristicEquationCheck, FrobeniusCharacteristicEquationExhaustiveReport,
    verify_frobenius_characteristic_equation_at_point,
    verify_frobenius_characteristic_equation_exhaustive,
};
pub use characteristic_polynomial::FrobeniusCharacteristicPolynomial;
pub use curve_type::{FrobeniusCurveType, FrobeniusCurveTypeReport};
pub use discriminant::FrobeniusDiscriminant;
pub use extension_counts::{
    FrobeniusExtensionCountReport, FrobeniusExtensionCountSequenceReport,
    FrobeniusExtensionEnumerationComparisonReport, compare_extension_count_with_enumeration,
};
pub use group_order::{GroupOrderReport, GroupOrderStrategy, HasseGroupOrderStrategy};
pub use hasse::{HasseBoundReport, verify_hasse_bound};
pub(crate) use hasse_multiple_search::hasse_multiple_search_report;
pub use hasse_multiple_search::{HasseMultipleSearchReport, HasseMultipleSearchStep};
pub use interval::HasseInterval;
pub use isogeny::{
    IsogenyFrobeniusRelation, IsogenyGraphFrobeniusReport, IsogenyGraphNodeFrobeniusData,
    verify_isogeny_frobenius_relation, verify_isogeny_graph_frobenius_relation,
};
pub use orbit::{
    FrobeniusOrbit, absolute_frobenius_orbit, absolute_frobenius_orbits_on_points,
    relative_frobenius_orbit, relative_frobenius_orbits_on_points,
};
pub use quadratic_twist::QuadraticTwistFrobeniusRelation;
pub use short_weierstrass::{
    absolute_frobenius_power_point, frobenius_twist_power, relative_frobenius_point,
};
pub use torsion::{
    FrobeniusOnExactTorsionPoint, FrobeniusOnExactTorsionReport,
    absolute_frobenius_on_exact_torsion, relative_frobenius_on_exact_torsion,
};
pub use torsion_matrix::{
    FrobeniusTorsionMatrixError, FrobeniusTorsionMatrixReport, ModNMatrix2, NTorsionBasis,
    frobenius_matrix_on_n_torsion_basis,
};
pub use trace::FrobeniusTrace;
pub use types::{AbsoluteFrobenius, RelativeFrobenius};
pub use zeta::FrobeniusLocalZetaFunction;

#[cfg(test)]
mod tests;
