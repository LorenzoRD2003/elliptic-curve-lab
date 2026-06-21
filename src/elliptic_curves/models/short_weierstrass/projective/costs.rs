use crate::elliptic_curves::CoordinateOperationCost;

/// Educational operation kinds for the current short-Weierstrass projective
/// baseline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ShortWeierstrassProjectiveOperationKind {
    FromAffine,
    ToAffine,
    Normalize,
    Neg,
    Add,
    Double,
    MixedAdd,
    ScalarMul,
}

/// Educational cost summary for one projective-side operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShortWeierstrassProjectiveOperationCost {
    kind: ShortWeierstrassProjectiveOperationKind,
    representation_cost: CoordinateOperationCost,
    affine_additions: usize,
    affine_doublings: usize,
    note: &'static str,
}

impl ShortWeierstrassProjectiveOperationCost {
    /// Builds one educational cost report for the current projective baseline.
    pub(crate) const fn new(
        kind: ShortWeierstrassProjectiveOperationKind,
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
    pub const fn kind(self) -> ShortWeierstrassProjectiveOperationKind {
        self.kind
    }

    /// Returns the counted coordinate work outside the reused affine group law.
    pub const fn representation_cost(self) -> CoordinateOperationCost {
        self.representation_cost
    }

    /// Returns how many affine additions are reused by this baseline route.
    pub const fn affine_additions(self) -> usize {
        self.affine_additions
    }

    /// Returns how many affine doublings are reused by this baseline route.
    pub const fn affine_doublings(self) -> usize {
        self.affine_doublings
    }

    /// Returns the educational note attached to the cost model.
    pub const fn note(self) -> &'static str {
        self.note
    }

    /// Returns the current baseline cost for one operation whose count does not
    /// depend on a caller-supplied scalar.
    pub const fn for_kind(kind: ShortWeierstrassProjectiveOperationKind) -> Self {
        match kind {
            ShortWeierstrassProjectiveOperationKind::FromAffine => from_affine_cost(),
            ShortWeierstrassProjectiveOperationKind::ToAffine => to_affine_cost(),
            ShortWeierstrassProjectiveOperationKind::Normalize => normalize_projective_cost(),
            ShortWeierstrassProjectiveOperationKind::Neg => neg_projective_cost(),
            ShortWeierstrassProjectiveOperationKind::Add => add_projective_cost(),
            ShortWeierstrassProjectiveOperationKind::Double => double_projective_cost(),
            ShortWeierstrassProjectiveOperationKind::MixedAdd => mixed_add_projective_cost(),
            ShortWeierstrassProjectiveOperationKind::ScalarMul => {
                ShortWeierstrassProjectiveOperationCost::new(
                    ShortWeierstrassProjectiveOperationKind::ScalarMul,
                    CoordinateOperationCost::new(0, 0, 0, 0),
                    0,
                    0,
                    "scalar multiplication now reuses the native projective add-and-double surface directly",
                )
            }
        }
    }

    /// Returns the current baseline scalar-multiplication cost model for one
    /// concrete non-negative scalar.
    pub fn for_scalar_mul(scalar: u64) -> Self {
        scalar_mul_projective_cost(scalar)
    }
}

/// Counts the current coordinate work needed to recover affine coordinates or to
/// normalize one finite representative in the public homogeneous chart
/// `x = X / Z`, `y = Y / Z`.
const fn normalization_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(0, 2, 0, 1)
}

/// Counts the cost of the public homogeneous-to-internal-Jacobian change of
/// variables
///
/// `X_J = X Z`, `Y_J = Y Z^2`, `Z_J = Z`.
const fn homogeneous_to_jacobian_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(0, 2, 1, 0)
}

/// Counts the cost of the internal Jacobian-to-public homogeneous change of
/// variables
///
/// `X = X_J Z_J`, `Y = Y_J`, `Z = Z_J^3`.
const fn jacobian_to_homogeneous_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(0, 2, 1, 0)
}

/// Counts the cost of the current Jacobian doubling core
///
/// `S = 4 X Y^2`, `M = 3 X^2 + a Z^4`,
///
/// `X3 = M^2 - 2S`, `Y3 = M(S - X3) - 8 Y^4`, `Z3 = 2 Y Z`.
const fn jacobian_double_core_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(12, 2, 6, 0)
}

/// Counts the cost of the current Jacobian addition core
///
/// `U1 = X1 Z2^2`, `U2 = X2 Z1^2`,
///
/// `S1 = Y1 Z2^3`, `S2 = Y2 Z1^3`,
///
/// `H = U2 - U1`, `I = (2H)^2`, `J = H I`,
///
/// `r = 2(S2 - S1)`, `V = U1 I`,
///
/// `X3 = r^2 - J - 2V`,
///
/// `Y3 = r(V - X3) - 2 S1 J`,
///
/// `Z3 = ((Z1 + Z2)^2 - Z1^2 - Z2^2) H`.
const fn jacobian_add_core_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(14, 8, 3, 0)
}

/// Counts the cost of the current Jacobian-affine mixed-add core
///
/// `U_2 = x_2 Z_1^2`, `S_2 = y_2 Z_1^3`,
///
/// `H = U_2 - X_1`, `I = (2H)^2`, `J = H I`,
///
/// `r = 2(S_2 - Y_1)`, `V = X_1 I`,
///
/// `X_3 = r^2 - J - 2V`,
///
/// `Y_3 = r(V - X_3) - 2 Y_1 J`,
///
/// `Z_3 = (Z_1 + H)^2 - Z_1^2 - H^2`.
const fn jacobian_mixed_add_core_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(9, 11, 5, 0)
}

/// Returns the current baseline cost for lifting a point from affine form.
const fn from_affine_cost() -> ShortWeierstrassProjectiveOperationCost {
    ShortWeierstrassProjectiveOperationCost::new(
        ShortWeierstrassProjectiveOperationKind::FromAffine,
        CoordinateOperationCost::new(0, 0, 0, 0),
        0,
        0,
        "the normalized chart reuses the affine coordinates directly and sets Z = 1",
    )
}

/// Returns the current baseline cost for recovering affine coordinates.
const fn to_affine_cost() -> ShortWeierstrassProjectiveOperationCost {
    ShortWeierstrassProjectiveOperationCost::new(
        ShortWeierstrassProjectiveOperationKind::ToAffine,
        normalization_cost(),
        0,
        0,
        "the current homogeneous chart divides both coordinates by the shared Z",
    )
}

/// Returns the current baseline cost for normalizing a finite point.
const fn normalize_projective_cost() -> ShortWeierstrassProjectiveOperationCost {
    ShortWeierstrassProjectiveOperationCost::new(
        ShortWeierstrassProjectiveOperationKind::Normalize,
        normalization_cost(),
        0,
        0,
        "normalization rescales one finite representative to the Z = 1 chart",
    )
}

/// Returns the current baseline cost for projective negation.
const fn neg_projective_cost() -> ShortWeierstrassProjectiveOperationCost {
    ShortWeierstrassProjectiveOperationCost::new(
        ShortWeierstrassProjectiveOperationKind::Neg,
        CoordinateOperationCost::new(0, 0, 0, 0),
        0,
        0,
        "the current counter tracks only additions, multiplications, squarings, and inversions",
    )
}

/// Returns the current baseline cost for adding two finite projective points
/// through the internal Jacobian engine.
const fn add_projective_cost() -> ShortWeierstrassProjectiveOperationCost {
    ShortWeierstrassProjectiveOperationCost::new(
        ShortWeierstrassProjectiveOperationKind::Add,
        homogeneous_to_jacobian_cost()
            .repeat(2)
            .combine(jacobian_add_core_cost())
            .combine(jacobian_to_homogeneous_cost()),
        0,
        0,
        "cost = two homogeneous-to-Jacobian chart changes + one Jacobian addition core + one Jacobian-to-homogeneous chart change",
    )
}

/// Returns the current baseline cost for doubling one finite projective point
/// through the internal Jacobian engine.
const fn double_projective_cost() -> ShortWeierstrassProjectiveOperationCost {
    ShortWeierstrassProjectiveOperationCost::new(
        ShortWeierstrassProjectiveOperationKind::Double,
        homogeneous_to_jacobian_cost()
            .combine(jacobian_double_core_cost())
            .combine(jacobian_to_homogeneous_cost()),
        0,
        0,
        "cost = one homogeneous-to-Jacobian chart change + one Jacobian doubling core + one Jacobian-to-homogeneous chart change",
    )
}

/// Returns the current baseline cost for mixing one projective and one affine
/// point.
const fn mixed_add_projective_cost() -> ShortWeierstrassProjectiveOperationCost {
    ShortWeierstrassProjectiveOperationCost::new(
        ShortWeierstrassProjectiveOperationKind::MixedAdd,
        homogeneous_to_jacobian_cost()
            .combine(jacobian_mixed_add_core_cost())
            .combine(jacobian_to_homogeneous_cost()),
        0,
        0,
        "cost = one homogeneous-to-Jacobian chart change + one Jacobian-affine mixed-add core + one Jacobian-to-homogeneous chart change",
    )
}

/// Returns the current baseline cost model for scalar multiplication.
fn scalar_mul_projective_cost(scalar: u64) -> ShortWeierstrassProjectiveOperationCost {
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

    ShortWeierstrassProjectiveOperationCost::new(
        ShortWeierstrassProjectiveOperationKind::ScalarMul,
        representation_cost,
        0,
        0,
        "scalar multiplication now counts repeated native projective additions and doublings",
    )
}
