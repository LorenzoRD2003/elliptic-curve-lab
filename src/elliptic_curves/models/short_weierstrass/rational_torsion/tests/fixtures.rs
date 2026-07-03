use num_bigint::BigInt;
use num_rational::BigRational;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::rational_torsion::{RationalTorsionGroup, RationalTorsionGroupShape},
};
use crate::fields::Q;

pub(super) struct RationalTorsionFixture {
    pub(super) name: &'static str,
    pub(super) curve: ShortWeierstrassCurve<Q>,
    pub(super) expected_group: RationalTorsionGroup,
    pub(super) sample_points: Vec<AffinePoint<Q>>,
}

pub(super) fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn point(x: i64, y: i64) -> AffinePoint<Q> {
    AffinePoint::new(q(x, 1), q(y, 1))
}

fn group(shape: RationalTorsionGroupShape) -> RationalTorsionGroup {
    RationalTorsionGroup::new(shape).expect("fixture torsion shape should satisfy Mazur")
}

/// `y² = x³ - x` has rational torsion `ℤ/2ℤ × ℤ/2ℤ`.
pub(super) fn product_two_two_fixture() -> RationalTorsionFixture {
    RationalTorsionFixture {
        name: "product-two-two",
        curve: ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("valid fixture curve"),
        expected_group: group(RationalTorsionGroupShape::ProductZ2Z2m { m: 1 }),
        sample_points: vec![
            AffinePoint::infinity(),
            point(-1, 0),
            point(0, 0),
            point(1, 0),
        ],
    }
}

/// `y² = x³ + 1` has rational torsion `ℤ/6ℤ`; `(2,3)` is a generator.
pub(super) fn cyclic_six_fixture() -> RationalTorsionFixture {
    RationalTorsionFixture {
        name: "cyclic-six",
        curve: ShortWeierstrassCurve::<Q>::new(q(0, 1), q(1, 1)).expect("valid fixture curve"),
        expected_group: group(RationalTorsionGroupShape::Cyclic { order: 6 }),
        sample_points: vec![
            AffinePoint::infinity(),
            point(-1, 0),
            point(0, 1),
            point(2, 3),
        ],
    }
}

/// `y² = x³ + x + 1` has trivial rational torsion.
pub(super) fn trivial_torsion_fixture() -> RationalTorsionFixture {
    RationalTorsionFixture {
        name: "trivial",
        curve: ShortWeierstrassCurve::<Q>::new(q(1, 1), q(1, 1)).expect("valid fixture curve"),
        expected_group: group(RationalTorsionGroupShape::Trivial),
        sample_points: vec![AffinePoint::infinity()],
    }
}

/// `y² = x³ - x/16` scales by `u = 2` to `y² = x³ - x`.
pub(super) fn rational_scaled_fixture() -> RationalTorsionFixture {
    RationalTorsionFixture {
        name: "scaled-product-two-two",
        curve: ShortWeierstrassCurve::<Q>::new(q(-1, 16), q(0, 1))
            .expect("valid rational fixture curve"),
        expected_group: group(RationalTorsionGroupShape::ProductZ2Z2m { m: 1 }),
        sample_points: vec![
            AffinePoint::infinity(),
            AffinePoint::new(q(-1, 4), q(0, 1)),
            AffinePoint::new(q(0, 1), q(0, 1)),
            AffinePoint::new(q(1, 4), q(0, 1)),
        ],
    }
}
