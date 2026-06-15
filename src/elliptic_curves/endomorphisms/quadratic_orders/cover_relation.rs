use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::quadratic_orders::ImaginaryQuadraticOrder;

/// One labeled cover relation in the Hasse diagram of candidate quadratic orders.
///
/// The orientation is from the larger order to the immediately contained
/// smaller one. The edge label is the relative index `[overorder : suborder]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuadraticOrderCoverRelation {
    overorder: ImaginaryQuadraticOrder,
    suborder: ImaginaryQuadraticOrder,
    index: BigUint,
}

impl QuadraticOrderCoverRelation {
    pub(crate) fn new(
        overorder: ImaginaryQuadraticOrder,
        suborder: ImaginaryQuadraticOrder,
        index: BigUint,
    ) -> Self {
        Self {
            overorder,
            suborder,
            index,
        }
    }

    /// Returns the larger order in the cover relation.
    pub fn overorder(&self) -> &ImaginaryQuadraticOrder {
        &self.overorder
    }

    /// Returns the immediately contained smaller order in the cover relation.
    pub fn suborder(&self) -> &ImaginaryQuadraticOrder {
        &self.suborder
    }

    /// Returns the relative index `[overorder : suborder]`.
    pub fn index(&self) -> &BigUint {
        &self.index
    }
}
