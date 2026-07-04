use num_rational::BigRational;

use crate::elliptic_curves::short_weierstrass::{
    group_law_core::{ShortWeierstrassFormulaPoint, ShortWeierstrassFormulaRunner},
    rational_torsion::{
        integral_model::{RationalIntegralModel, integral_rational_to_bigint},
        mazur::MAZUR_CYCLIC_ORDERS,
        reduction_mod_p::{
            formula_ops::ReductionFormulaOps,
            good_prime::GoodReductionPrime,
            small_prime_field::{ReductionPrime, ReductionResidue},
        },
    },
};

/// Short-Weierstrass point on the reduced curve over `ℤ/pℤ`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum ReducedPoint {
    Infinity,
    Finite {
        x: ReductionResidue,
        y: ReductionResidue,
    },
}

/// Reduced point whose order is compatible with Mazur's theorem over `Q`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct MazurReducedPoint {
    point: ReducedPoint,
    order: usize,
}

/// Integral short-Weierstrass model reduced modulo a good prime.
///
/// It stores `E_p: y² = x³ + (A mod p)x + (B mod p)` and the
/// good-prime certificate `Δ(E_int) ≠ 0 mod p`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct ReducedShortWeierstrassCurve {
    good_prime: GoodReductionPrime,
    a: ReductionResidue,
    b: ReductionResidue,
}

impl ReducedShortWeierstrassCurve {
    /// Reduces an integral short-Weierstrass model modulo the first good prime.
    ///
    /// Complexity: linear in the number of primes tested by the good-prime
    /// search, plus two coefficient reductions.
    pub(super) fn from_integral_model(model: &RationalIntegralModel) -> Option<Self> {
        let good_prime = GoodReductionPrime::first_for_integral_model(model)?;
        let prime = good_prime.prime();
        let curve = model.curve();
        let a = reduce_integral_coefficient(prime, curve.a(), "A");
        let b = reduce_integral_coefficient(prime, curve.b(), "B");

        Some(Self { good_prime, a, b })
    }

    pub(super) fn good_prime(&self) -> GoodReductionPrime {
        self.good_prime
    }

    pub(super) fn prime(&self) -> ReductionPrime {
        self.good_prime.prime()
    }

    pub(super) fn a(&self) -> ReductionResidue {
        self.a
    }

    pub(super) fn b(&self) -> ReductionResidue {
        self.b
    }

    /// Enumerates all points of `E_p(𝔽_p)` by scanning every affine pair.
    ///
    /// Complexity: `Θ(p²)` modular operations.
    pub(super) fn points(&self) -> Vec<ReducedPoint> {
        let prime = self.prime();
        let mut points = vec![ReducedPoint::Infinity];

        for x in prime.residues() {
            let rhs = self.rhs(x);
            for y in prime.residues() {
                if prime.mul(y, y) == rhs {
                    points.push(ReducedPoint::Finite { x, y });
                }
            }
        }
        points
    }

    pub(super) fn contains(&self, point: ReducedPoint) -> bool {
        match point {
            ReducedPoint::Infinity => true,
            ReducedPoint::Finite { x, y } => {
                let prime = self.prime();
                prime.mul(y, y) == self.rhs(x)
            }
        }
    }

    pub(super) fn neg(&self, point: ReducedPoint) -> ReducedPoint {
        match point {
            ReducedPoint::Infinity => ReducedPoint::Infinity,
            ReducedPoint::Finite { x, y } => ReducedPoint::Finite {
                x,
                y: self.prime().neg(y),
            },
        }
    }

    /// Adds two points on `E_p: y² = x³ + Ax + B`.
    ///
    /// Complexity: `Θ(1)` modular operations.
    pub(super) fn add(&self, left: ReducedPoint, right: ReducedPoint) -> ReducedPoint {
        let ops = ReductionFormulaOps::new(self.prime());
        let runner = ShortWeierstrassFormulaRunner::new(&ops, &self.a);
        ReducedPoint::from_formula_point(
            runner
                .add_points(&left.to_formula_point(), &right.to_formula_point())
                .expect("reduced affine formula addition should not fail"),
        )
    }

    /// Multiplies a reduced point by a non-negative scalar.
    ///
    /// Complexity: `Θ(log n)` reduced-curve additions.
    pub(super) fn mul_scalar(&self, point: ReducedPoint, n: usize) -> ReducedPoint {
        let ops = ReductionFormulaOps::new(self.prime());
        let runner = ShortWeierstrassFormulaRunner::new(&ops, &self.a);
        ReducedPoint::from_formula_point(
            runner
                .mul_point(&point.to_formula_point(), n)
                .expect("reduced affine formula scalar multiplication should not fail"),
        )
    }

    /// Returns the first Mazur-permitted order killing `point`, if any.
    ///
    /// Complexity: constant-time with respect to the curve, since Mazur's list
    /// has fixed size; each scalar multiplication costs `Θ(log m)`.
    pub(super) fn mazur_order(&self, point: ReducedPoint) -> Option<usize> {
        if !self.contains(point) {
            return None;
        }
        if point == ReducedPoint::Infinity {
            return Some(1);
        }

        for order in MAZUR_CYCLIC_ORDERS {
            if self.mul_scalar(point, *order) == ReducedPoint::Infinity {
                return Some(*order);
            }
        }
        None
    }

    /// Filters the reduced point set to points killed by a Mazur-permitted
    /// order.
    ///
    /// Complexity: `Θ(N)`, where `N = #E_p(𝔽_p)`, because the Mazur order list
    /// is fixed.
    pub(super) fn mazur_order_points(&self) -> Vec<MazurReducedPoint> {
        self.points()
            .into_iter()
            .filter_map(|point| {
                self.mazur_order(point)
                    .map(|order| MazurReducedPoint { point, order })
            })
            .collect()
    }

    fn rhs(&self, x: ReductionResidue) -> ReductionResidue {
        let prime = self.prime();
        let x_sq = prime.mul(x, x);
        let x_cube = prime.mul(x_sq, x);
        let ax = prime.mul(self.a, x);
        prime.add(prime.add(x_cube, ax), self.b)
    }
}

impl ReducedPoint {
    fn to_formula_point(self) -> ShortWeierstrassFormulaPoint<ReductionResidue> {
        match self {
            Self::Infinity => ShortWeierstrassFormulaPoint::Infinity,
            Self::Finite { x, y } => ShortWeierstrassFormulaPoint::Affine { x, y },
        }
    }

    fn from_formula_point(point: ShortWeierstrassFormulaPoint<ReductionResidue>) -> Self {
        match point {
            ShortWeierstrassFormulaPoint::Infinity => Self::Infinity,
            ShortWeierstrassFormulaPoint::Affine { x, y } => Self::Finite { x, y },
        }
    }
}

impl MazurReducedPoint {
    pub(super) fn point(&self) -> ReducedPoint {
        self.point
    }

    pub(super) fn order(&self) -> usize {
        self.order
    }
}

fn reduce_integral_coefficient(
    prime: ReductionPrime,
    coefficient: &BigRational,
    name: &'static str,
) -> ReductionResidue {
    let coefficient = integral_rational_to_bigint(coefficient).unwrap_or_else(|_| {
        panic!("RationalIntegralModel should have integral coefficient {name}")
    });
    prime.reduce_bigint(&coefficient)
}
