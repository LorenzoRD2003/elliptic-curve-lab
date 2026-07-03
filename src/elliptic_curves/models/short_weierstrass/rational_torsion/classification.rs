use crate::elliptic_curves::short_weierstrass::rational_torsion::mazur::{
    MAZUR_CYCLIC_ORDERS, MAZUR_PRODUCT_PARAMETERS,
};

use super::RationalTorsionError;

/// Candidate Mazur-shape classification for `E(Q)_tors`.
///
/// Mazur's theorem says that a torsion subgroup over `Q` is either trivial,
/// nontrivially cyclic `ℤ/nℤ`, where `2 ≤ n ≤ 10` or `n = 12`, or a product
/// `ℤ/2ℤ × ℤ/2mℤ`, where `1 ≤ m ≤ 4`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RationalTorsionGroupShape {
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
            Self::Cyclic { order } if MAZUR_CYCLIC_ORDERS.contains(&order) => Ok(()),
            Self::Cyclic { order } => Err(RationalTorsionError::InvalidMazurShape {
                family: "cyclic",
                value: order,
            }),
            Self::ProductZ2Z2m { m } if MAZUR_PRODUCT_PARAMETERS.contains(&m) => Ok(()),
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
pub struct RationalTorsionGroup {
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

    /// Classifies a verified finite rational-torsion set from exact point orders.
    ///
    /// The input is expected to contain the exact order of every verified point,
    /// including `1` for `O`. If the resulting finite group does not match a
    /// Mazur shape, this returns an internal consistency error.
    pub(crate) fn from_verified_point_orders(
        point_orders: &[usize],
    ) -> Result<Self, RationalTorsionError> {
        let point_count = point_orders.len();
        let max_order = point_orders.iter().copied().max().unwrap_or(0);

        let shape = if point_count == 1 {
            RationalTorsionGroupShape::Trivial
        } else if max_order == point_count {
            RationalTorsionGroupShape::Cyclic { order: point_count }
        } else if point_count.is_multiple_of(4) && max_order * 2 == point_count {
            RationalTorsionGroupShape::ProductZ2Z2m { m: point_count / 4 }
        } else {
            return Err(RationalTorsionError::InconsistentMazurShape { point_count });
        };

        Self::new(shape).map_err(|_| RationalTorsionError::InconsistentMazurShape { point_count })
    }

    /// Returns the stored Mazur shape.
    pub fn shape(self) -> RationalTorsionGroupShape {
        self.shape
    }

    /// Returns the cardinality of the classified torsion group.
    pub fn cardinality(self) -> usize {
        match self.shape {
            RationalTorsionGroupShape::Trivial => 1,
            RationalTorsionGroupShape::Cyclic { order } => order,
            RationalTorsionGroupShape::ProductZ2Z2m { m } => 4 * m,
        }
    }
}
