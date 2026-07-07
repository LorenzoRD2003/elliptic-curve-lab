use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::{ShortestFloorPathReport, volcano_structure::VolcanoStructureRole},
};

/// Certified structural data for one stored graph node.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VolcanoStructureNodeReport {
    node_id: IsogenyGraphNodeId,
    distance_to_floor: usize,
    level: usize,
    role: VolcanoStructureRole,
    floor_path: ShortestFloorPathReport,
}

impl VolcanoStructureNodeReport {
    pub(crate) fn from_floor_path(
        certified_depth: usize,
        floor_path: ShortestFloorPathReport,
    ) -> Self {
        let distance_to_floor = floor_path.distance_to_floor();
        let level = certified_depth
            .checked_sub(distance_to_floor)
            .expect("certified depth is the maximum certified distance");
        let role = VolcanoStructureRole::for_level(level, certified_depth);

        Self {
            node_id: floor_path.start(),
            distance_to_floor,
            level,
            role,
            floor_path,
        }
    }

    /// Returns the stored graph node.
    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.node_id
    }

    /// Returns the certified distance `δ(v) = dist(v, V_d)`.
    pub fn distance_to_floor(&self) -> usize {
        self.distance_to_floor
    }

    /// Returns the certified level `d̂ - δ(v)`.
    pub fn level(&self) -> usize {
        self.level
    }

    /// Returns this node's structural role.
    pub fn role(&self) -> VolcanoStructureRole {
        self.role
    }

    /// Returns the shortest floor path certifying `δ(v)`.
    pub fn floor_path(&self) -> &ShortestFloorPathReport {
        &self.floor_path
    }
}
