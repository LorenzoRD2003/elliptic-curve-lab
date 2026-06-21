use crate::elliptic_curves::CoordinateOperationCost;

/// Educational operation kinds for the current general-Weierstrass projective
/// baseline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeneralWeierstrassProjectiveOperationKind {
    FromAffine,
    ToAffine,
    Normalize,
    Neg,
    Add,
    Double,
    MixedAdd,
    ScalarMul,
}

/// Educational cost summary for one general-Weierstrass projective operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GeneralWeierstrassProjectiveOperationCost {
    kind: GeneralWeierstrassProjectiveOperationKind,
    representation_cost: CoordinateOperationCost,
    affine_additions: usize,
    affine_doublings: usize,
    note: &'static str,
}

impl GeneralWeierstrassProjectiveOperationCost {
    pub(crate) const fn new(
        kind: GeneralWeierstrassProjectiveOperationKind,
        representation_cost: CoordinateOperationCost,
        affine_additions: usize,
        affine_doublings: usize,
        note: &'static str,
    ) -> Self {
        Self {
            kind,
            representation_cost,
            affine_additions,
            affine_doublings,
            note,
        }
    }

    /// Returns the operation kind.
    pub const fn kind(self) -> GeneralWeierstrassProjectiveOperationKind {
        self.kind
    }

    /// Returns the counted coordinate work outside the reused affine group law.
    pub const fn representation_cost(self) -> CoordinateOperationCost {
        self.representation_cost
    }

    /// Returns how many affine additions are reused by this staged route.
    pub const fn affine_additions(self) -> usize {
        self.affine_additions
    }

    /// Returns how many affine doublings are reused by this staged route.
    pub const fn affine_doublings(self) -> usize {
        self.affine_doublings
    }

    /// Returns the educational note attached to the cost model.
    pub const fn note(self) -> &'static str {
        self.note
    }

    /// Returns the current baseline cost for one operation whose count does not
    /// depend on a caller-supplied scalar.
    pub const fn for_kind(kind: GeneralWeierstrassProjectiveOperationKind) -> Self {
        match kind {
            GeneralWeierstrassProjectiveOperationKind::FromAffine => from_affine_cost(),
            GeneralWeierstrassProjectiveOperationKind::ToAffine => to_affine_cost(),
            GeneralWeierstrassProjectiveOperationKind::Normalize => normalize_projective_cost(),
            GeneralWeierstrassProjectiveOperationKind::Neg => neg_projective_cost(),
            GeneralWeierstrassProjectiveOperationKind::Add => add_projective_cost(),
            GeneralWeierstrassProjectiveOperationKind::Double => double_projective_cost(),
            GeneralWeierstrassProjectiveOperationKind::MixedAdd => mixed_add_projective_cost(),
            GeneralWeierstrassProjectiveOperationKind::ScalarMul => Self::new(
                GeneralWeierstrassProjectiveOperationKind::ScalarMul,
                CoordinateOperationCost::new(0, 0, 0, 0),
                0,
                0,
                "scalar multiplication now reuses the native projective add-and-double surface directly",
            ),
        }
    }

    /// Returns the current baseline scalar-multiplication cost model for one
    /// concrete non-negative scalar.
    pub fn for_scalar_mul(scalar: u64) -> Self {
        scalar_mul_projective_cost(scalar)
    }
}

const fn normalization_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(0, 2, 0, 1)
}

/// Counts the cost of preparing the homogenized secant slope data
///
/// `Y2 Z1 - Y1 Z2`, `X2 Z1 - X1 Z2`,
///
/// together with the shared denominator data `Z1 Z2`, `Z1^2`, `Z2^2`,
/// `(X2 Z1 - X1 Z2)^2`.
const fn general_add_slope_setup_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(2, 5, 3, 0)
}

/// Counts the cost of assembling the homogenized `X3` numerator for the
/// current addition formula.
const fn general_add_x_recovery_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(8, 6, 1, 0)
}

/// Counts the cost of assembling the homogenized `Y3` numerator and the final
/// shared projective denominator for the current addition formula.
const fn general_add_y_recovery_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(10, 7, 0, 0)
}

/// Counts the cost of preparing the homogenized secant slope data for the
/// mixed case `Z_2 = 1`
///
/// `y_2 Z_1 - Y_1`, `x_2 Z_1 - X_1`,
///
/// together with `Z_1^2` and `(x_2 Z_1 - X_1)^2`.
const fn general_mixed_add_slope_setup_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(2, 2, 2, 0)
}

/// Counts the cost of assembling the homogenized `X3` numerator for the
/// current mixed-add formula.
const fn general_mixed_add_x_recovery_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(4, 8, 1, 0)
}

/// Counts the cost of assembling the homogenized `Y3` numerator and the final
/// shared projective denominator for the current mixed-add formula.
const fn general_mixed_add_y_recovery_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(4, 10, 1, 0)
}

/// Counts the cost of preparing the homogenized tangent slope data
///
/// `N = 3 X^2 + 2 a2 X Z + a4 Z^2 - a1 Y Z`,
///
/// `D = 2 Y + a1 X + a3 Z`,
///
/// together with `Z^2` and `D^2`.
const fn general_double_tangent_setup_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(4, 7, 2, 0)
}

/// Counts the cost of assembling the homogenized `X3` numerator for the
/// current doubling formula.
const fn general_double_x_recovery_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(5, 5, 1, 0)
}

/// Counts the cost of assembling the homogenized `Y3` numerator and the final
/// shared projective denominator for the current doubling formula.
const fn general_double_y_recovery_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(6, 3, 1, 0)
}

const fn from_affine_cost() -> GeneralWeierstrassProjectiveOperationCost {
    GeneralWeierstrassProjectiveOperationCost::new(
        GeneralWeierstrassProjectiveOperationKind::FromAffine,
        CoordinateOperationCost::new(0, 0, 0, 0),
        0,
        0,
        "the normalized chart reuses the affine coordinates directly and sets Z = 1",
    )
}

const fn to_affine_cost() -> GeneralWeierstrassProjectiveOperationCost {
    GeneralWeierstrassProjectiveOperationCost::new(
        GeneralWeierstrassProjectiveOperationKind::ToAffine,
        normalization_cost(),
        0,
        0,
        "the current homogeneous chart divides both coordinates by the shared Z",
    )
}

const fn normalize_projective_cost() -> GeneralWeierstrassProjectiveOperationCost {
    GeneralWeierstrassProjectiveOperationCost::new(
        GeneralWeierstrassProjectiveOperationKind::Normalize,
        normalization_cost(),
        0,
        0,
        "normalization rescales one finite representative to the Z = 1 chart",
    )
}

const fn neg_projective_cost() -> GeneralWeierstrassProjectiveOperationCost {
    GeneralWeierstrassProjectiveOperationCost::new(
        GeneralWeierstrassProjectiveOperationKind::Neg,
        CoordinateOperationCost::new(0, 2, 0, 0),
        0,
        0,
        "projective negation updates Y by the model-specific involution while keeping X and Z",
    )
}

const fn add_projective_cost() -> GeneralWeierstrassProjectiveOperationCost {
    GeneralWeierstrassProjectiveOperationCost::new(
        GeneralWeierstrassProjectiveOperationKind::Add,
        general_add_slope_setup_cost()
            .combine(general_add_x_recovery_cost())
            .combine(general_add_y_recovery_cost()),
        0,
        0,
        "cost = homogenized secant-slope setup + homogenized X3 assembly + homogenized Y3 and shared-denominator assembly",
    )
}

const fn double_projective_cost() -> GeneralWeierstrassProjectiveOperationCost {
    GeneralWeierstrassProjectiveOperationCost::new(
        GeneralWeierstrassProjectiveOperationKind::Double,
        general_double_tangent_setup_cost()
            .combine(general_double_x_recovery_cost())
            .combine(general_double_y_recovery_cost()),
        0,
        0,
        "cost = homogenized tangent-slope setup + homogenized X3 assembly + homogenized Y3 and shared-denominator assembly",
    )
}

const fn mixed_add_projective_cost() -> GeneralWeierstrassProjectiveOperationCost {
    GeneralWeierstrassProjectiveOperationCost::new(
        GeneralWeierstrassProjectiveOperationKind::MixedAdd,
        general_mixed_add_slope_setup_cost()
            .combine(general_mixed_add_x_recovery_cost())
            .combine(general_mixed_add_y_recovery_cost()),
        0,
        0,
        "cost = homogenized mixed secant-slope setup with Z_2 = 1 + homogenized X3 assembly + homogenized Y3 and shared-denominator assembly",
    )
}

fn scalar_mul_projective_cost(scalar: u64) -> GeneralWeierstrassProjectiveOperationCost {
    let bit_length = u64::BITS as usize - scalar.leading_zeros() as usize;
    let projective_doublings = bit_length.saturating_sub(1);
    let projective_additions = scalar.count_ones() as usize;
    let representation_cost = add_projective_cost()
        .representation_cost()
        .repeat(projective_additions)
        .combine(
            double_projective_cost()
                .representation_cost()
                .repeat(projective_doublings),
        );

    GeneralWeierstrassProjectiveOperationCost::new(
        GeneralWeierstrassProjectiveOperationKind::ScalarMul,
        representation_cost,
        0,
        0,
        "scalar multiplication now counts repeated native projective additions and doublings",
    )
}
