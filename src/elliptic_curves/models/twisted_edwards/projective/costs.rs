use crate::elliptic_curves::CoordinateOperationCost;

/// Educational operation kinds for the current twisted-Edwards projective
/// engine in extended `(X:Y:Z:T)` coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TwistedEdwardsProjectiveOperationKind {
    FromAffine,
    ToAffine,
    Neg,
    Add,
    Double,
    MixedAdd,
    ScalarMul,
}

/// Educational cost summary for one twisted-Edwards projective operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TwistedEdwardsProjectiveOperationCost {
    kind: TwistedEdwardsProjectiveOperationKind,
    representation_cost: CoordinateOperationCost,
    affine_additions: usize,
    affine_doublings: usize,
    note: &'static str,
}

impl TwistedEdwardsProjectiveOperationCost {
    pub(crate) const fn new(
        kind: TwistedEdwardsProjectiveOperationKind,
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
    pub const fn kind(self) -> TwistedEdwardsProjectiveOperationKind {
        self.kind
    }

    /// Returns the counted coordinate work outside any delegated affine route.
    pub const fn representation_cost(self) -> CoordinateOperationCost {
        self.representation_cost
    }

    /// Returns how many affine additions are delegated by this route.
    pub const fn affine_additions(self) -> usize {
        self.affine_additions
    }

    /// Returns how many affine doublings are delegated by this route.
    pub const fn affine_doublings(self) -> usize {
        self.affine_doublings
    }

    /// Returns the educational note attached to this cost model.
    pub const fn note(self) -> &'static str {
        self.note
    }

    /// Returns the current baseline cost for one operation whose count does not
    /// depend on a caller-supplied scalar.
    pub const fn for_kind(kind: TwistedEdwardsProjectiveOperationKind) -> Self {
        match kind {
            TwistedEdwardsProjectiveOperationKind::FromAffine => from_affine_cost(),
            TwistedEdwardsProjectiveOperationKind::ToAffine => to_affine_cost(),
            TwistedEdwardsProjectiveOperationKind::Neg => neg_projective_cost(),
            TwistedEdwardsProjectiveOperationKind::Add => add_projective_cost(),
            TwistedEdwardsProjectiveOperationKind::Double => double_projective_cost(),
            TwistedEdwardsProjectiveOperationKind::MixedAdd => mixed_add_projective_cost(),
            TwistedEdwardsProjectiveOperationKind::ScalarMul => Self::new(
                TwistedEdwardsProjectiveOperationKind::ScalarMul,
                CoordinateOperationCost::new(0, 0, 0, 0),
                0,
                0,
                "scalar multiplication reuses the native projective add-and-double surface directly",
            ),
        }
    }

    /// Returns the current baseline scalar-multiplication cost model for one
    /// concrete non-negative scalar.
    pub fn for_scalar_mul(scalar: u64) -> Self {
        scalar_mul_projective_cost(scalar)
    }
}

const fn from_affine_cost() -> TwistedEdwardsProjectiveOperationCost {
    TwistedEdwardsProjectiveOperationCost::new(
        TwistedEdwardsProjectiveOperationKind::FromAffine,
        CoordinateOperationCost::new(0, 1, 0, 0),
        0,
        0,
        "embedding into the normalized extended chart sets Z = 1 and computes T = x y",
    )
}

const fn to_affine_cost() -> TwistedEdwardsProjectiveOperationCost {
    TwistedEdwardsProjectiveOperationCost::new(
        TwistedEdwardsProjectiveOperationKind::ToAffine,
        CoordinateOperationCost::new(0, 2, 0, 1),
        0,
        0,
        "affine recovery divides X and Y by the shared Z coordinate",
    )
}

const fn neg_projective_cost() -> TwistedEdwardsProjectiveOperationCost {
    TwistedEdwardsProjectiveOperationCost::new(
        TwistedEdwardsProjectiveOperationKind::Neg,
        CoordinateOperationCost::new(0, 0, 0, 0),
        0,
        0,
        "the current counter tracks only additions, multiplications, squarings, and inversions; sign flips are not charged separately",
    )
}

const fn extended_add_core_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(7, 9, 0, 0)
}

const fn extended_double_core_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(6, 6, 3, 0)
}

const fn extended_mixed_add_core_cost() -> CoordinateOperationCost {
    CoordinateOperationCost::new(6, 8, 0, 0)
}

const fn add_projective_cost() -> TwistedEdwardsProjectiveOperationCost {
    TwistedEdwardsProjectiveOperationCost::new(
        TwistedEdwardsProjectiveOperationKind::Add,
        extended_add_core_cost(),
        0,
        0,
        "cost = native extended addition in (X:Y:Z:T) coordinates with no affine delegation",
    )
}

const fn double_projective_cost() -> TwistedEdwardsProjectiveOperationCost {
    TwistedEdwardsProjectiveOperationCost::new(
        TwistedEdwardsProjectiveOperationKind::Double,
        extended_double_core_cost(),
        0,
        0,
        "cost = native extended doubling in (X:Y:Z:T) coordinates with no affine delegation",
    )
}

const fn mixed_add_projective_cost() -> TwistedEdwardsProjectiveOperationCost {
    TwistedEdwardsProjectiveOperationCost::new(
        TwistedEdwardsProjectiveOperationKind::MixedAdd,
        extended_mixed_add_core_cost(),
        0,
        0,
        "cost = native mixed addition specialized to an affine right input already in the Z_2 = 1 chart",
    )
}

fn scalar_mul_projective_cost(scalar: u64) -> TwistedEdwardsProjectiveOperationCost {
    if scalar == 0 {
        return TwistedEdwardsProjectiveOperationCost::new(
            TwistedEdwardsProjectiveOperationKind::ScalarMul,
            CoordinateOperationCost::new(0, 0, 0, 0),
            0,
            0,
            "the zero scalar returns the projective identity immediately",
        );
    }

    let bits = u64::BITS as usize - scalar.leading_zeros() as usize;
    let additions = scalar.count_ones() as usize;
    let doublings = bits.saturating_sub(1);

    TwistedEdwardsProjectiveOperationCost::new(
        TwistedEdwardsProjectiveOperationKind::ScalarMul,
        add_projective_cost()
            .representation_cost()
            .repeat(additions)
            .combine(
                double_projective_cost()
                    .representation_cost()
                    .repeat(doublings),
            ),
        0,
        0,
        "cost = one native projective addition per set bit plus one native projective doubling per remaining processed bit",
    )
}
