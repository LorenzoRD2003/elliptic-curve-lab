use super::RationalTorsionError;

/// Candidate Mazur-shape classification for `E(Q)_tors`.
///
/// Mazur's theorem says that a torsion subgroup over `Q` is either cyclic
/// `ℤ/nℤ`, where `1 ≤ n ≤ 10` or `n = 12`, or a product
/// `ℤ/2ℤ × ℤ/2mℤ`, where `1 ≤ m ≤ 4`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RationalTorsionGroupShape {
    /// The trivial group `{O}`.
    Trivial,
    /// A cyclic group `ℤ/nℤ`.
    Cyclic {
        /// The order `n` of the cyclic torsion group.
        order: usize,
    },
    /// A product group `ℤ/2ℤ × ℤ/2mℤ`.
    ProductZ2Z2m {
        /// The integer `m` in the Mazur product family.
        m: usize,
    },
}

impl RationalTorsionGroupShape {
    fn validate(self) -> Result<(), RationalTorsionError> {
        match self {
            Self::Trivial => Ok(()),
            Self::Cyclic { order } if matches!(order, 1..=10 | 12) => Ok(()),
            Self::Cyclic { order } => Err(RationalTorsionError::InvalidMazurShape {
                family: "cyclic",
                value: order,
            }),
            Self::ProductZ2Z2m { m } if matches!(m, 1..=4) => Ok(()),
            Self::ProductZ2Z2m { m } => Err(RationalTorsionError::InvalidMazurShape {
                family: "product",
                value: m,
            }),
        }
    }
}

/// Validated Mazur-shape classification for the rational torsion group
/// `E(Q)_tors`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RationalTorsionGroup {
    shape: RationalTorsionGroupShape,
}

impl RationalTorsionGroup {
    /// Builds a validated rational torsion group classification.
    ///
    /// This is the only constructor because `RationalTorsionGroup` represents
    /// an already-certified Mazur shape, not an arbitrary candidate.
    pub(crate) fn new(shape: RationalTorsionGroupShape) -> Result<Self, RationalTorsionError> {
        shape.validate()?;
        Ok(Self { shape })
    }

    /// Returns the stored Mazur shape.
    pub(crate) fn shape(self) -> RationalTorsionGroupShape {
        self.shape
    }

    /// Returns the cardinality of the classified torsion group.
    pub(crate) fn cardinality(self) -> usize {
        match self.shape {
            RationalTorsionGroupShape::Trivial => 1,
            RationalTorsionGroupShape::Cyclic { order } => order,
            RationalTorsionGroupShape::ProductZ2Z2m { m } => 4 * m,
        }
    }
}
