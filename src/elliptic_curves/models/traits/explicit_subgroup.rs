use std::hash::Hash;

use crate::elliptic_curves::{CurveError, traits::GroupCurveModel};

/// Crate-private capability for exhaustively validating an explicit finite
/// subgroup of curve points.
///
/// This lives with curve-model traits rather than under `isogenies`, because
/// the mathematical story is group-theoretic before it is isogeny-specific:
/// we are checking that a set of points forms a subgroup of the ambient
/// additive curve group.
pub(crate) trait ExplicitSubgroupCurveModel: GroupCurveModel
where
    Self::Point: Clone + Eq + Hash,
{
    fn validate_explicit_point_subgroup(
        &self,
        points: &HashSet<Self::Point>,
    ) -> Result<(), CurveError> {
        if points.is_empty() {
            return Err(CurveError::GroupAxiomViolation {
                axiom: "explicit subgroup must be non-empty",
            });
        }

        let identity = self.identity();
        if !points.contains(&identity) {
            return Err(CurveError::GroupAxiomViolation {
                axiom: "explicit subgroup must contain the identity",
            });
        }

        if points.iter().any(|point| !self.contains(point)) {
            return Err(CurveError::PointNotOnCurve);
        }

        for point in points {
            let inverse = self.neg(point);
            if !points.contains(&inverse) {
                return Err(CurveError::GroupAxiomViolation {
                    axiom: "explicit subgroup must be closed under negation",
                });
            }
        }

        for left in points {
            for right in points {
                let sum = self.add(left, right)?;
                if !points.contains(&sum) {
                    return Err(CurveError::GroupAxiomViolation {
                        axiom: "explicit subgroup must be closed under addition",
                    });
                }
            }
        }

        Ok(())
    }
}

impl<C: GroupCurveModel> ExplicitSubgroupCurveModel for C where C::Point: Clone + Eq + Hash {}
use std::collections::HashSet;
