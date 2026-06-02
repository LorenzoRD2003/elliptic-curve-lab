use std::hash::Hash;

use super::{
    EdgeTargetWitness, IsogenyGraph, IsogenyGraphEdge, IsogenyGraphEdgeId, IsogenyGraphError,
};
use crate::elliptic_curves::{
    CurveIsomorphism, ShortWeierstrassCurve, ShortWeierstrassIsomorphism,
};
use crate::fields::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    ComposedIsogeny, Isogeny, IsogenyError, IsomorphismIsogeny, ScalarMultiplicationIsogeny,
    VeluIsogeny, VerifiableIsogeny, maps_equal_exhaustively,
};

type Curve<F> = ShortWeierstrassCurve<F>;
type TargetTransport<F> = IsomorphismIsogeny<ShortWeierstrassIsomorphism<F>>;
type EdgeMap<F> =
    ComposedIsogeny<VeluIsogeny<Curve<F>>, TargetTransport<F>, Curve<F>, Curve<F>, Curve<F>>;

/// Status of the reverse-direction graph edge corresponding to one directed edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReverseEdgeStatus {
    Missing,
    PresentButNotVerifiedAsDual,
    VerifiedAsDual,
}

/// Exhaustive small-field local verification summary for an educational isogeny graph.
#[derive(Clone, Debug)]
pub struct IsogenyGraphVerificationReport {
    pub checked_edges: usize,
    pub edges_mapping_domain_to_codomain: usize,
    pub edges_mapping_kernel_to_identity: usize,
    pub edges_homomorphism_ok: usize,
    pub reverse_edges_verified_as_dual: usize,
    pub reverse_edge_statuses: Vec<(IsogenyGraphEdgeId, ReverseEdgeStatus)>,
}

impl<F> IsogenyGraph<ShortWeierstrassCurve<F>>
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    /// Exhaustively verifies each stored edge as a small rational-point map.
    ///
    /// This verification is intentionally local and graph-aware:
    ///
    /// - it reconstructs the stored edge as `transport ∘ Vélu`
    /// - it checks that the map lands on the declared target representative
    /// - it checks that kernel points map to the codomain identity
    /// - it checks the additive homomorphism law on rational points
    /// - it inspects reverse-direction edges already present in the graph and
    ///   classifies them as missing, present-but-unverified, or verified duals
    ///
    /// It does not run a fresh global dual search for every edge.
    pub fn verify_locally(&self) -> Result<IsogenyGraphVerificationReport, IsogenyGraphError> {
        let mut report = IsogenyGraphVerificationReport {
            checked_edges: 0,
            edges_mapping_domain_to_codomain: 0,
            edges_mapping_kernel_to_identity: 0,
            edges_homomorphism_ok: 0,
            reverse_edges_verified_as_dual: 0,
            reverse_edge_statuses: Vec::with_capacity(self.edge_count()),
        };

        for edge in self.edges() {
            report.checked_edges += 1;

            let phi = self.reconstruct_edge_map(edge)?;

            match phi.verify_maps_domain_to_codomain() {
                Ok(()) => report.edges_mapping_domain_to_codomain += 1,
                Err(IsogenyError::ImagePointNotOnCodomain) => {}
                Err(other) => return Err(other.into()),
            }

            match phi.verify_maps_kernel_to_identity() {
                Ok(()) => report.edges_mapping_kernel_to_identity += 1,
                Err(IsogenyError::KernelPointDoesNotMapToIdentity) => {}
                Err(other) => return Err(other.into()),
            }

            match phi.verify_homomorphism() {
                Ok(()) => report.edges_homomorphism_ok += 1,
                Err(IsogenyError::HomomorphismViolation) => {}
                Err(other) => return Err(other.into()),
            }

            let reverse_status = self.reverse_edge_status(edge)?;
            if reverse_status == ReverseEdgeStatus::VerifiedAsDual {
                report.reverse_edges_verified_as_dual += 1;
            }
            report
                .reverse_edge_statuses
                .push((edge.id(), reverse_status));
        }

        Ok(report)
    }

    fn reconstruct_edge_map(
        &self,
        edge: &IsogenyGraphEdge<ShortWeierstrassCurve<F>>,
    ) -> Result<EdgeMap<F>, IsogenyGraphError> {
        let source_curve = self
            .node(edge.source())
            .ok_or(IsogenyGraphError::MissingSourceNode(edge.source()))?
            .representative()
            .clone();
        let target_curve = self
            .node(edge.target())
            .ok_or(IsogenyGraphError::MissingTargetNode(edge.target()))?
            .representative()
            .clone();
        let generator = edge
            .kernel()
            .points()
            .get(1)
            .ok_or(IsogenyGraphError::InvalidDegree)?
            .clone();
        let velu = VeluIsogeny::from_generator(source_curve, generator)?;

        if velu.kernel() != edge.kernel() {
            return Err(IsogenyError::KernelMismatch.into());
        }

        let transport_isomorphism =
            self.edge_target_isomorphism(edge, velu.codomain(), &target_curve)?;
        let transport = IsomorphismIsogeny::new(transport_isomorphism);

        ComposedIsogeny::new_strict(velu, transport).map_err(Into::into)
    }

    fn edge_target_isomorphism(
        &self,
        edge: &IsogenyGraphEdge<ShortWeierstrassCurve<F>>,
        raw_codomain: &ShortWeierstrassCurve<F>,
        target_curve: &ShortWeierstrassCurve<F>,
    ) -> Result<ShortWeierstrassIsomorphism<F>, IsogenyGraphError> {
        let isomorphism = match edge.target_witness() {
            EdgeTargetWitness::Identity => {
                ShortWeierstrassIsomorphism::new(raw_codomain.clone(), F::one())
                    .map_err(IsogenyError::from)
                    .map_err(IsogenyGraphError::from)?
            }
            EdgeTargetWitness::Explicit(isomorphism) => isomorphism.clone(),
        };

        if isomorphism.domain() != raw_codomain || isomorphism.codomain() != target_curve {
            return Err(IsogenyError::CompositionDomainCodomainMismatch.into());
        }

        Ok(isomorphism)
    }

    fn reverse_edge_status(
        &self,
        edge: &IsogenyGraphEdge<ShortWeierstrassCurve<F>>,
    ) -> Result<ReverseEdgeStatus, IsogenyGraphError> {
        let reverse_candidates = self
            .outgoing_edges(edge.target())
            .filter(|candidate| candidate.target() == edge.source())
            .filter(|candidate| candidate.degree() == edge.degree())
            .collect::<Vec<_>>();

        if reverse_candidates.is_empty() {
            return Ok(ReverseEdgeStatus::Missing);
        }

        for reverse_edge in reverse_candidates {
            if self.reverse_edge_is_verified_dual(edge, reverse_edge)? {
                return Ok(ReverseEdgeStatus::VerifiedAsDual);
            }
        }

        Ok(ReverseEdgeStatus::PresentButNotVerifiedAsDual)
    }

    fn reverse_edge_is_verified_dual(
        &self,
        edge: &IsogenyGraphEdge<ShortWeierstrassCurve<F>>,
        reverse_edge: &IsogenyGraphEdge<ShortWeierstrassCurve<F>>,
    ) -> Result<bool, IsogenyGraphError> {
        let phi = self.reconstruct_edge_map(edge)?;
        let psi = self.reconstruct_edge_map(reverse_edge)?;
        let degree = u64::try_from(edge.degree()).expect("tiny educational degrees should fit");

        let left_composition =
            ComposedIsogeny::new_strict(phi, psi).map_err(IsogenyGraphError::from)?;
        let left_scalar = ScalarMultiplicationIsogeny::new(
            self.node(edge.source())
                .ok_or(IsogenyGraphError::MissingSourceNode(edge.source()))?
                .representative()
                .clone(),
            degree,
        )?;
        if !maps_equal_exhaustively::<_, _, Curve<F>, Curve<F>>(&left_composition, &left_scalar)? {
            return Ok(false);
        }

        let psi = self.reconstruct_edge_map(reverse_edge)?;
        let phi = self.reconstruct_edge_map(edge)?;
        let right_composition =
            ComposedIsogeny::new_strict(psi, phi).map_err(IsogenyGraphError::from)?;
        let right_scalar = ScalarMultiplicationIsogeny::new(
            self.node(edge.target())
                .ok_or(IsogenyGraphError::MissingTargetNode(edge.target()))?
                .representative()
                .clone(),
            degree,
        )?;

        maps_equal_exhaustively::<_, _, Curve<F>, Curve<F>>(&right_composition, &right_scalar)
            .map_err(Into::into)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::{Field, Fp};
    use crate::isogenies::graphs::IsogenyGraphBuilder;
    use crate::isogenies::{Isogeny, VerifiableIsogeny};

    type F17 = Fp<17>;
    type Curve17 = ShortWeierstrassCurve<F17>;

    fn f17_graph() -> crate::isogenies::graphs::IsogenyGraph<Curve17> {
        let curve = Curve17::new(F17::from_i64(1), F17::from_i64(0)).expect("valid curve");
        IsogenyGraphBuilder::new(curve, 2)
            .max_depth(3)
            .build()
            .expect("F17 graph should build")
    }

    #[test]
    fn every_graph_edge_reconstructs_a_valid_velu_isogeny() {
        let graph = f17_graph();

        for edge in graph.edges() {
            let reconstructed = graph
                .reconstruct_edge_map(edge)
                .expect("stored graph edge should reconstruct");

            assert_eq!(reconstructed.degree(), edge.degree());
            assert!(reconstructed.verify_maps_domain_to_codomain().is_ok());
            assert!(reconstructed.verify_maps_kernel_to_identity().is_ok());
            assert!(reconstructed.verify_homomorphism().is_ok());
        }
    }
}
