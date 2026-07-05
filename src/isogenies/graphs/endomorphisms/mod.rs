mod edge_relation;
mod graph_endomorphism_report;
mod node_candidates;
pub mod refinement;

#[cfg(test)]
mod volcano_report;

pub use edge_relation::{IsogenyEdgeEndomorphismReport, IsogenyEdgeEndomorphismTentativeRelation};
pub use graph_endomorphism_report::{
    IsogenyGraphEndomorphismEdgeReport, IsogenyGraphEndomorphismNodeReport,
    IsogenyGraphEndomorphismReport,
};

#[cfg(test)]
pub(crate) use volcano_report::{EndomorphismVolcanoReport, VolcanoHeuristicComparison};
