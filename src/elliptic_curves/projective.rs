use core::fmt;
use std::hash::{Hash, Hasher};

use crate::elliptic_curves::{AffinePoint, CurveError};
use crate::fields::traits::Field;

/// Educational projective point representation with an explicit infinity variant.
///
/// Finite points are stored as triples `(X : Y : Z)` with `Z != 0`, interpreted
/// through the homogeneous affine recovery `x = X / Z`, `y = Y / Z`.
///
/// The distinguished point at infinity is modeled explicitly instead of being
/// encoded by a degenerate `Z = 0` triple.
pub enum ProjectivePoint<F: Field> {
    Infinity,
    Finite { x: F::Elem, y: F::Elem, z: F::Elem },
}

/// Small educational count of coordinate operations.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CoordinateOperationCost {
    additions: usize,
    multiplications: usize,
    squarings: usize,
    inversions: usize,
}

impl CoordinateOperationCost {
    /// Builds one educational operation-cost value object.
    pub const fn new(
        additions: usize,
        multiplications: usize,
        squarings: usize,
        inversions: usize,
    ) -> Self {
        Self {
            additions,
            multiplications,
            squarings,
            inversions,
        }
    }

    /// Returns the counted additions.
    pub const fn additions(self) -> usize {
        self.additions
    }

    /// Returns the counted multiplications.
    pub const fn multiplications(self) -> usize {
        self.multiplications
    }

    /// Returns the counted squarings.
    pub const fn squarings(self) -> usize {
        self.squarings
    }

    /// Returns the counted inversions.
    pub const fn inversions(self) -> usize {
        self.inversions
    }

    /// Returns the pointwise sum of two operation counts.
    pub const fn combine(self, other: Self) -> Self {
        Self::new(
            self.additions + other.additions,
            self.multiplications + other.multiplications,
            self.squarings + other.squarings,
            self.inversions + other.inversions,
        )
    }

    /// Scales all counters by the same non-negative integer.
    pub const fn repeat(self, count: usize) -> Self {
        Self::new(
            self.additions * count,
            self.multiplications * count,
            self.squarings * count,
            self.inversions * count,
        )
    }

    /// Returns whether every counter is zero.
    pub const fn is_zero(self) -> bool {
        self.additions == 0
            && self.multiplications == 0
            && self.squarings == 0
            && self.inversions == 0
    }
}

impl<F: Field> ProjectivePoint<F> {
    /// Builds a finite projective point without validating curve membership.
    pub fn new(x: F::Elem, y: F::Elem, z: F::Elem) -> Self {
        Self::Finite { x, y, z }
    }

    /// Builds the distinguished projective point at infinity.
    pub fn infinity() -> Self {
        Self::Infinity
    }

    /// Lifts an affine point into the normalized projective chart `Z = 1`.
    pub fn from_affine(point: &AffinePoint<F>) -> Self {
        match point {
            AffinePoint::Infinity => Self::Infinity,
            AffinePoint::Finite { x, y } => Self::Finite {
                x: x.clone(),
                y: y.clone(),
                z: F::one(),
            },
        }
    }

    /// Returns whether this point is the distinguished identity element.
    pub fn is_identity(&self) -> bool {
        matches!(self, Self::Infinity)
    }

    /// Returns whether the stored finite representative already lies in the
    /// normalized `Z = 1` chart.
    ///
    /// The point at infinity is treated as already normalized.
    pub fn is_normalized(&self) -> bool {
        match self {
            Self::Infinity => true,
            Self::Finite { z, .. } => F::eq(z, &F::one()),
        }
    }

    /// Recovers the corresponding affine point.
    ///
    /// Finite points require `Z` to be invertible; a stored `Z = 0` therefore
    /// returns the underlying field inversion failure.
    pub fn to_affine(&self) -> Result<AffinePoint<F>, CurveError> {
        match self {
            Self::Infinity => Ok(AffinePoint::Infinity),
            Self::Finite { x, y, z } => {
                let z_inverse = F::inverse(z).map_err(CurveError::Field)?;
                Ok(AffinePoint::new(
                    F::mul(x, &z_inverse),
                    F::mul(y, &z_inverse),
                ))
            }
        }
    }

    /// Returns an equivalent projective representative with `Z = 1`.
    pub fn normalize(&self) -> Result<Self, CurveError> {
        match self {
            Self::Infinity => Ok(Self::Infinity),
            Self::Finite { x, y, z } => {
                let z_inverse = F::inverse(z).map_err(CurveError::Field)?;
                Ok(Self::Finite {
                    x: F::mul(x, &z_inverse),
                    y: F::mul(y, &z_inverse),
                    z: F::one(),
                })
            }
        }
    }

    /// Transports the stored coordinates through a caller-supplied map.
    pub fn map_coordinates<G, M>(&self, mut map: M) -> ProjectivePoint<G>
    where
        G: Field,
        M: FnMut(&F::Elem) -> G::Elem,
    {
        match self {
            Self::Infinity => ProjectivePoint::Infinity,
            Self::Finite { x, y, z } => ProjectivePoint::new(map(x), map(y), map(z)),
        }
    }

    /// Returns a compact coordinate string or `O` for the identity.
    pub fn to_coordinates_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        match self {
            Self::Infinity => "O".to_string(),
            Self::Finite { x, y, z } => format!("({x} : {y} : {z})"),
        }
    }

    /// Returns whether both values store exactly the same representative.
    ///
    /// This is stricter than `PartialEq`: projectively equivalent finite
    /// representatives such as `(X : Y : Z)` and `(λX : λY : λZ)` compare equal
    /// under `==`, but they do not have the same stored representative unless
    /// their coordinates also match componentwise.
    pub fn has_same_representative_as(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Infinity, Self::Infinity) => true,
            (
                Self::Finite {
                    x: x_left,
                    y: y_left,
                    z: z_left,
                },
                Self::Finite {
                    x: x_right,
                    y: y_right,
                    z: z_right,
                },
            ) => F::eq(x_left, x_right) && F::eq(y_left, y_right) && F::eq(z_left, z_right),
            _ => false,
        }
    }
}

impl<F: Field> PartialEq for ProjectivePoint<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Infinity, Self::Infinity) => true,
            (
                Self::Finite {
                    x: x_left,
                    y: y_left,
                    z: z_left,
                },
                Self::Finite {
                    x: x_right,
                    y: y_right,
                    z: z_right,
                },
            ) => {
                if F::is_zero(z_left) || F::is_zero(z_right) {
                    return self.has_same_representative_as(other);
                }

                F::eq(&F::mul(x_left, z_right), &F::mul(x_right, z_left))
                    && F::eq(&F::mul(y_left, z_right), &F::mul(y_right, z_left))
            }
            _ => false,
        }
    }
}

impl<F: Field> Eq for ProjectivePoint<F> {}

impl<F: Field> Hash for ProjectivePoint<F>
where
    F::Elem: Hash,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            Self::Infinity => {
                0_u8.hash(state);
            }
            Self::Finite { x, y, z } => {
                if F::is_zero(z) {
                    1_u8.hash(state);
                    x.hash(state);
                    y.hash(state);
                    z.hash(state);
                    return;
                }

                2_u8.hash(state);
                match self.normalize() {
                    Ok(Self::Finite { x, y, .. }) => {
                        x.hash(state);
                        y.hash(state);
                    }
                    Ok(Self::Infinity) => {
                        0_u8.hash(state);
                    }
                    Err(_) => {
                        x.hash(state);
                        y.hash(state);
                        z.hash(state);
                    }
                }
            }
        }
    }
}

impl<F: Field> Clone for ProjectivePoint<F> {
    fn clone(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Finite { x, y, z } => Self::Finite {
                x: x.clone(),
                y: y.clone(),
                z: z.clone(),
            },
        }
    }
}

impl<F: Field> fmt::Display for ProjectivePoint<F>
where
    F::Elem: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_coordinates_string())
    }
}

impl<F: Field> fmt::Debug for ProjectivePoint<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Infinity => write!(f, "ProjectivePoint::Infinity"),
            Self::Finite { x, y, z } => f
                .debug_struct("ProjectivePoint")
                .field("x", x)
                .field("y", y)
                .field("z", z)
                .finish(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    use crate::elliptic_curves::{AffinePoint, CoordinateOperationCost, ProjectivePoint};
    use crate::fields::{Fp, traits::Field};

    type F7 = Fp<7>;

    #[test]
    fn affine_points_lift_into_the_normalized_projective_chart() {
        let affine = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5));
        let projective = ProjectivePoint::from_affine(&affine);

        assert_eq!(
            projective,
            ProjectivePoint::new(F7::from_i64(2), F7::from_i64(5), F7::one())
        );
        assert!(projective.is_normalized());
    }

    #[test]
    fn projective_to_affine_divides_by_the_shared_z_coordinate() {
        let point = ProjectivePoint::<F7>::new(F7::from_i64(6), F7::from_i64(1), F7::from_i64(3));

        assert_eq!(
            point.to_affine(),
            Ok(AffinePoint::new(F7::from_i64(2), F7::from_i64(5)))
        );
    }

    #[test]
    fn normalization_rescales_finite_points_to_z_equal_one() {
        let point = ProjectivePoint::<F7>::new(F7::from_i64(6), F7::from_i64(1), F7::from_i64(3));
        let normalized = point.normalize().expect("nonzero z should normalize");

        assert_eq!(
            normalized,
            ProjectivePoint::new(F7::from_i64(2), F7::from_i64(5), F7::one())
        );
        assert!(normalized.is_normalized());
    }

    #[test]
    fn equality_uses_projective_equivalence_not_raw_representative_equality() {
        let left = ProjectivePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5), F7::one());
        let right = ProjectivePoint::<F7>::new(F7::from_i64(6), F7::from_i64(1), F7::from_i64(3));

        assert_eq!(left, right);
        assert!(!left.has_same_representative_as(&right));
    }

    #[test]
    fn projectively_equal_finite_points_hash_the_same_way() {
        let left = ProjectivePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5), F7::one());
        let right = ProjectivePoint::<F7>::new(F7::from_i64(6), F7::from_i64(1), F7::from_i64(3));

        let mut left_hasher = DefaultHasher::new();
        left.hash(&mut left_hasher);

        let mut right_hasher = DefaultHasher::new();
        right.hash(&mut right_hasher);

        assert_eq!(left_hasher.finish(), right_hasher.finish());
    }

    #[test]
    fn finite_points_with_zero_z_fail_affine_recovery_and_normalization() {
        let point = ProjectivePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5), F7::zero());
        let same_raw = ProjectivePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5), F7::zero());
        let different_raw =
            ProjectivePoint::<F7>::new(F7::from_i64(4), F7::from_i64(3), F7::zero());

        assert!(point.to_affine().is_err());
        assert!(point.normalize().is_err());
        assert!(!point.is_normalized());
        assert_eq!(point, same_raw);
        assert_ne!(point, different_raw);
    }

    #[test]
    fn map_coordinates_transforms_all_finite_components() {
        let point = ProjectivePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5), F7::from_i64(3));
        let image = point.map_coordinates::<F7, _>(|coordinate| F7::add(coordinate, &F7::one()));

        assert_eq!(
            image,
            ProjectivePoint::new(F7::from_i64(3), F7::from_i64(6), F7::from_i64(4))
        );
    }

    #[test]
    fn coordinate_cost_combine_and_repeat_accumulate_counters_componentwise() {
        let left = CoordinateOperationCost {
            additions: 1,
            multiplications: 2,
            squarings: 3,
            inversions: 4,
        };
        let right = CoordinateOperationCost {
            additions: 3,
            multiplications: 4,
            squarings: 5,
            inversions: 6,
        };

        assert_eq!(
            left.combine(right),
            CoordinateOperationCost::new(4, 6, 8, 10)
        );
        assert_eq!(left.repeat(2), CoordinateOperationCost::new(2, 4, 6, 8));
        assert!(CoordinateOperationCost::default().is_zero());
    }
}
