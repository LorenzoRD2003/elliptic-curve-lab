use crate::elliptic_curves::endomorphisms::{
    ImaginaryQuadraticOrder, ImaginaryQuadraticOrderError, QuadraticDiscriminant,
    QuadraticDiscriminantFactorization, QuadraticOrderIndexError,
};
use crate::numerics::positive_divisors;
use num_bigint::BigUint;

/// Candidate imaginary quadratic orders satisfying `ℤ[π] ⊆ O_f ⊆ O_K`.
///
/// If the Frobenius discriminant factors as `Δ_π = v^2 D_K`, then every
/// intermediate quadratic order in the same imaginary quadratic field has
/// conductor `f` dividing `v`. This value object enumerates exactly those
/// candidate conductors and packages the corresponding orders.
///
/// Mathematically, these candidates carry the divisibility poset on the
/// conductors `f | v`, equivalently the order-containment poset on the
/// corresponding `O_f`. That poset is not a total chain in general:
/// for example, when `v = 6`, the conductors `2` and `3` are both valid but
/// incomparable. The current representation is intentionally smaller: it stores
/// the candidates as one list sorted by increasing conductor, while leaving any
/// future Hasse-diagram or poset view as a derived presentation layer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EndomorphismRingCandidateSet {
    factorization: QuadraticDiscriminantFactorization,
    candidate_orders: Vec<ImaginaryQuadraticOrder>,
}

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

impl QuadraticDiscriminantFactorization {
    /// Enumerates the candidate orders `O_f` with `ℤ[π] ⊆ O_f ⊆ O_K`.
    ///
    /// If `Δ_π = v^2 D_K`, then the possible conductors are exactly the
    /// positive divisors `f | v`.
    ///
    /// The returned orders are sorted by increasing conductor, so the first
    /// element is the maximal order `O_K` and the last element is the
    /// Frobenius-generated order `ℤ[π]`.
    ///
    /// This ordering is only a deterministic storage convention. It should **not**
    /// be read as saying that the candidate orders form a total chain under
    /// inclusion for every conductor `v`.
    ///
    /// Complexity: dominated by `num-prime`.
    pub fn endomorphism_ring_candidates(
        &self,
    ) -> Result<EndomorphismRingCandidateSet, ImaginaryQuadraticOrderError> {
        EndomorphismRingCandidateSet::from_factorization(self.clone())
    }
}

impl EndomorphismRingCandidateSet {
    /// Builds the candidate-order set from the canonical factorization
    /// `Δ_π = v^2 D_K`.
    ///
    /// Complexity: dominated by `num-prime`, plus `Θ(τ(v))` order constructions.
    pub fn from_factorization(
        factorization: QuadraticDiscriminantFactorization,
    ) -> Result<Self, ImaginaryQuadraticOrderError> {
        let fundamental_discriminant = factorization.fundamental_discriminant().clone();
        let candidate_orders = positive_divisors(factorization.conductor())
            .into_iter()
            .map(|conductor| {
                ImaginaryQuadraticOrder::new(fundamental_discriminant.clone(), conductor)
            })
            .collect::<Result<Vec<_>, _>>()?;

        debug_assert!(!candidate_orders.is_empty());

        Ok(Self {
            factorization,
            candidate_orders,
        })
    }

    /// Builds the candidate-order set from the discriminant `Δ_π`.
    ///
    /// Complexity: dominated by `num-prime`.
    pub fn from_discriminant(
        discriminant: &QuadraticDiscriminant,
    ) -> Result<Self, ImaginaryQuadraticOrderError> {
        let factorization = discriminant
            .factorization()
            .map_err(|_| ImaginaryQuadraticOrderError::NonImaginaryOrderDiscriminant)?;
        Self::from_factorization(factorization)
    }

    /// Returns the canonical factorization `Δ_π = v^2 D_K`.
    pub fn factorization(&self) -> &QuadraticDiscriminantFactorization {
        &self.factorization
    }

    /// Returns the original discriminant `Δ_π`.
    pub fn discriminant(&self) -> &QuadraticDiscriminant {
        self.factorization.discriminant()
    }

    /// Returns the fundamental discriminant `D_K`.
    pub fn fundamental_discriminant(&self) -> &QuadraticDiscriminant {
        self.factorization.fundamental_discriminant()
    }

    /// Returns the positive square root factor `v` from `Δ_π = v^2 D_K`.
    pub fn frobenius_conductor(&self) -> &BigUint {
        self.factorization.conductor()
    }

    /// Returns the candidate orders `O_f`, sorted by increasing conductor `f`.
    ///
    /// This is a deterministic list view of the candidates. The mathematically
    /// natural relation between them is the divisibility poset on conductors,
    /// equivalently the containment poset on orders.
    pub fn candidate_orders(&self) -> &[ImaginaryQuadraticOrder] {
        &self.candidate_orders
    }

    /// Returns the number of candidate orders.
    pub fn len(&self) -> usize {
        self.candidate_orders.len()
    }

    /// Returns whether the candidate list is empty.
    ///
    /// For a valid factorization this is always `false`, since `1 | v` for
    /// every positive conductor `v`.
    pub fn is_empty(&self) -> bool {
        self.candidate_orders.is_empty()
    }

    /// Returns the maximal order `O_K`, corresponding to `f = 1`.
    pub fn maximal_order(&self) -> &ImaginaryQuadraticOrder {
        self.candidate_orders
            .first()
            .expect("candidate-order sets are always non-empty")
    }

    /// Returns the Frobenius-generated order `ℤ[π]`, corresponding to `f = v`.
    pub fn frobenius_order(&self) -> &ImaginaryQuadraticOrder {
        self.candidate_orders
            .last()
            .expect("candidate-order sets are always non-empty")
    }

    /// Returns the relative index `[candidate_order : ℤ[π]]`.
    ///
    /// If `Δ_π = v^2 D_K`, the Frobenius-generated order is `ℤ[π] = O_v`.
    /// For any candidate `O_f` with `f | v`, this index is `v / f`.
    ///
    /// Complexity: `Θ(1)` big-integer arithmetic.
    pub fn index_over_frobenius_order(
        &self,
        candidate_order: &ImaginaryQuadraticOrder,
    ) -> Result<BigUint, QuadraticOrderIndexError> {
        candidate_order.index_of_suborder(self.frobenius_order())
    }

    /// Returns the labeled cover relations in the Hasse diagram of the candidate-order poset.
    ///
    /// The poset order is inclusion of quadratic orders, equivalently reverse
    /// divisibility on conductors. An edge `O_g -> O_f` appears exactly when
    /// `O_f ⊆ O_g` and there is no intermediate candidate order strictly
    /// between them in the current candidate set.
    ///
    /// Complexity: `Θ(n^3)` order comparisons in the current straightforward
    /// implementation, where `n = self.len()`.
    pub fn hasse_cover_relations(&self) -> Vec<QuadraticOrderCoverRelation> {
        let mut relations = Vec::new();

        for (over_index, overorder) in self.candidate_orders.iter().enumerate() {
            for (sub_index, suborder) in self.candidate_orders.iter().enumerate() {
                if over_index >= sub_index
                    || !suborder.is_suborder_of(overorder)
                    || overorder == suborder
                {
                    continue;
                }

                let has_intermediate =
                    self.candidate_orders
                        .iter()
                        .enumerate()
                        .any(|(middle_index, middle_order)| {
                            middle_index != over_index
                                && middle_index != sub_index
                                && suborder.is_suborder_of(middle_order)
                                && middle_order.is_suborder_of(overorder)
                                && middle_order != overorder
                                && middle_order != suborder
                        });

                if !has_intermediate {
                    relations.push(QuadraticOrderCoverRelation {
                        overorder: overorder.clone(),
                        suborder: suborder.clone(),
                        index: overorder
                            .index_of_suborder(suborder)
                            .expect("Hasse edges stay inside one candidate-order poset"),
                    });
                }
            }
        }

        relations
    }
}

impl QuadraticOrderCoverRelation {
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
