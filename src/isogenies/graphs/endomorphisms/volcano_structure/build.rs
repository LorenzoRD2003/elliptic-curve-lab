use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph,
    endomorphisms::{
        VolcanoSearchError,
        volcano_structure::{
            UncertifiedVolcanoNodeReport, VolcanoStructureReport, VolcanoStructureUncertifiedReason,
        },
    },
};
use crate::numerics::validate_positive_prime;

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    /// Builds a certified structural report for this stored `ℓ`-volcano graph.
    ///
    /// The method tries to certify `δ(v) = dist(v, V_d)` for every stored node
    /// using [`Self::find_shortest_floor_path`]. Certified nodes are grouped
    /// into levels by `level(v) = d̂ - δ(v)`, where `d̂` is the largest
    /// certified distance in the stored graph. Nodes whose distance cannot be
    /// certified remain in [`VolcanoStructureReport::uncertified_nodes`].
    ///
    /// Complexity: `Θ(Σ_v δ(v))` certified path steps over stored nodes, plus
    /// linear grouping by certified level. The method does not derive
    /// endomorphism rings and does not build new isogeny graphs.
    pub fn volcano_structure_report(
        &self,
        prime: &BigUint,
    ) -> Result<VolcanoStructureReport, VolcanoSearchError> {
        validate_positive_prime(prime)?;

        let mut floor_paths = Vec::new();
        let mut uncertified_nodes = Vec::new();

        for node in self.nodes() {
            match self.find_shortest_floor_path(node.id(), prime) {
                Ok(floor_path) => floor_paths.push(floor_path),
                Err(error) => {
                    let reason = VolcanoStructureUncertifiedReason::from_search_error(error)?;
                    uncertified_nodes.push(UncertifiedVolcanoNodeReport::new(node.id(), reason));
                }
            }
        }

        Ok(VolcanoStructureReport::from_floor_paths(
            prime.clone(),
            floor_paths,
            uncertified_nodes,
        ))
    }
}
