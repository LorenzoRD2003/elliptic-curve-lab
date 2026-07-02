use crate::fields::traits::Field;

/// Extension trait for field families that can decide and produce square roots.
///
/// This trait is intentionally small so later APIs can ask only for the extra
/// capability they need, rather than assuming every field backend has a usable
/// square-root algorithm.
///
/// Important scope note:
///
/// This trait does **not** promise that every implementation uses a complete
/// or asymptotically efficient algorithm. In the current educational stage of
/// the crate, some backends may use simple baseline procedures as long as the
/// behavior is documented honestly.
///
/// Contract:
///
/// - `sqrt(x)` returns `Some(y)` only when `y^2 = x`
/// - when square roots are not available in the field, it returns `None`
/// - if two square roots exist, the implementation may return either one
pub trait SqrtField: Field {
    /// Returns one square root of `x` when it exists in the field.
    fn sqrt(x: &Self::Elem) -> Option<Self::Elem>;

    /// Returns whether `x` admits a square root in the field.
    fn has_square_root(x: &Self::Elem) -> bool {
        Self::sqrt(x).is_some()
    }

    /// Returns one square root together with its additive inverse.
    ///
    /// In characteristic different from `2`, these are exactly the two roots
    /// of `y^2 = x`, except when `x = 0`, where both entries are zero.
    fn sqrt_pair(x: &Self::Elem) -> Option<(Self::Elem, Self::Elem)> {
        let root = Self::sqrt(x)?;
        Some((root.clone(), Self::neg(&root)))
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::traits::*;

    use crate::fields::traits::SqrtField;

    type F5 = crate::fields::Fp5;

    #[test]
    fn has_square_root_matches_quadratic_residues_over_f5() {
        assert!(F5::has_square_root(&F5::zero()));
        assert!(F5::has_square_root(&F5::one()));
        assert!(F5::has_square_root(&F5::from_i64(-1)));
        assert!(!F5::has_square_root(&F5::from_i64(2)));
    }

    #[test]
    fn sqrt_returns_a_value_whose_square_is_the_input() {
        let root = F5::sqrt(&F5::from_i64(4)).expect("4 should be a square in F5");

        assert!(F5::eq(&F5::square(&root), &F5::from_i64(4)));
    }

    #[test]
    fn sqrt_pair_returns_opposite_roots() {
        let (left, right) = F5::sqrt_pair(&F5::one()).expect("1 should be a square in F5");

        assert!(F5::eq(&F5::square(&left), &F5::one()));
        assert!(F5::eq(&F5::square(&right), &F5::one()));
        assert!(F5::eq(&right, &F5::neg(&left)));
    }
}
