use crate::fields::traits::Field;

/// Extension trait for field families that can decide and produce cube roots.
///
/// This trait parallels [`crate::fields::SqrtField`], but for the equation
/// `y^3 = x`.
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
/// - `cbrt(x)` returns `Some(y)` only when `y^3 = x`
/// - when cube roots are not available in the field, it returns `None`
/// - if multiple cube roots exist, the implementation may return any one of
///   them
pub trait CbrtField: Field {
    /// Returns one cube root of `x` when it exists in the field.
    fn cbrt(x: &Self::Elem) -> Option<Self::Elem>;

    /// Returns whether `x` admits a cube root in the field.
    fn has_cube_root(x: &Self::Elem) -> bool {
        Self::cbrt(x).is_some()
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::CbrtField;
    use crate::fields::{Field, Fp};

    type F7 = Fp<7>;

    #[test]
    fn has_cube_root_matches_cubic_residues_over_f7() {
        assert!(F7::has_cube_root(&F7::zero()));
        assert!(F7::has_cube_root(&F7::one()));
        assert!(F7::has_cube_root(&F7::from_i64(-1)));
        assert!(!F7::has_cube_root(&F7::from_i64(2)));
    }

    #[test]
    fn cbrt_returns_a_value_whose_cube_is_the_input() {
        let root = F7::cbrt(&F7::from_i64(-1)).expect("-1 should be a cube in F7");

        assert!(F7::eq(&F7::cube(&root), &F7::from_i64(-1)));
    }
}
