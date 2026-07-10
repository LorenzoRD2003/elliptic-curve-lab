use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, QuadraticClassGroup},
    quadratic_ideals::PrimeNormIdeal,
    quadratic_orders::QuadraticDiscriminant,
};

use super::{bu, maximal_order, z};

#[derive(Clone, Copy)]
struct GpIdealFormFixture {
    discriminant: i64,
    ell: u64,
    root: u64,
    raw: (i64, i64, i64),
    reduced: (i64, i64, i64),
    ramified: bool,
}

impl GpIdealFormFixture {
    fn ideal(self) -> PrimeNormIdeal {
        let order = maximal_order(self.discriminant);
        if self.ramified {
            PrimeNormIdeal::ramified(order, bu(self.ell))
                .expect("GP fixture should use a ramified prime ideal")
        } else {
            PrimeNormIdeal::split(order, bu(self.ell), bu(self.root))
                .expect("GP fixture should use a split prime ideal")
        }
    }

    fn raw_form(self) -> BinaryQuadraticForm {
        form(self.raw)
    }

    fn reduced_form(self) -> BinaryQuadraticForm {
        form(self.reduced)
    }
}

fn form((a, b, c): (i64, i64, i64)) -> BinaryQuadraticForm {
    BinaryQuadraticForm::new(z(a), z(b), z(c))
}

fn gp_ideal_form_fixtures() -> Vec<GpIdealFormFixture> {
    // External source: GP/PARI 2.17 via:
    //
    // ideal_form(D,ell,r)=my(m=2*ell,b);
    //   for(k=0,m-1,if((k-D)%2==0 && (k-r)%ell==0,b=k;break));
    //   my(c=(b^2-D)/(4*ell)); [D,ell,r,b,c,qfbred(Qfb(ell,b,c))]
    //
    // Here `r` is the root of `x² ≡ D (mod ℓ)` stored by `PrimeNormIdeal`.
    vec![
        GpIdealFormFixture {
            discriminant: -20,
            ell: 3,
            root: 1,
            raw: (3, 4, 3),
            reduced: (2, 2, 3),
            ramified: false,
        },
        GpIdealFormFixture {
            discriminant: -20,
            ell: 3,
            root: 2,
            raw: (3, 2, 2),
            reduced: (2, 2, 3),
            ramified: false,
        },
        GpIdealFormFixture {
            discriminant: -23,
            ell: 3,
            root: 1,
            raw: (3, 1, 2),
            reduced: (2, -1, 3),
            ramified: false,
        },
        GpIdealFormFixture {
            discriminant: -23,
            ell: 3,
            root: 2,
            raw: (3, 5, 4),
            reduced: (2, 1, 3),
            ramified: false,
        },
        GpIdealFormFixture {
            discriminant: -31,
            ell: 5,
            root: 2,
            raw: (5, 7, 4),
            reduced: (2, -1, 4),
            ramified: false,
        },
        GpIdealFormFixture {
            discriminant: -31,
            ell: 5,
            root: 3,
            raw: (5, 3, 2),
            reduced: (2, 1, 4),
            ramified: false,
        },
        GpIdealFormFixture {
            discriminant: -84,
            ell: 5,
            root: 1,
            raw: (5, 6, 6),
            reduced: (5, 4, 5),
            ramified: false,
        },
        GpIdealFormFixture {
            discriminant: -84,
            ell: 5,
            root: 4,
            raw: (5, 4, 5),
            reduced: (5, 4, 5),
            ramified: false,
        },
        GpIdealFormFixture {
            discriminant: -23,
            ell: 23,
            root: 0,
            raw: (23, 23, 6),
            reduced: (1, 1, 6),
            ramified: true,
        },
    ]
}

#[test]
fn gp_fixtures_for_prime_norm_ideal_form_labels_have_expected_reductions() {
    for fixture in gp_ideal_form_fixtures() {
        let ideal = fixture.ideal();
        let class_group =
            QuadraticClassGroup::new(QuadraticDiscriminant::new(fixture.discriminant))
                .expect("GP fixture discriminants should define class groups");
        let raw_form = fixture.raw_form();
        let reduced_form = fixture.reduced_form();

        assert_eq!(ideal.norm(), &bu(fixture.ell));
        assert_eq!(ideal.root_mod_ell(), &bu(fixture.root));
        assert_eq!(raw_form.a(), &z(fixture.ell as i64));
        assert_eq!(raw_form.discriminant(), z(fixture.discriminant));
        assert!(raw_form.is_primitive());
        assert!(raw_form.is_positive_definite());
        assert_eq!(
            raw_form
                .reduce_positive_definite()
                .expect("GP raw forms are positive definite"),
            reduced_form
        );
        assert_eq!(reduced_form.discriminant(), z(fixture.discriminant));
        assert!(reduced_form.is_reduced_positive_definite());
        assert!(
            class_group
                .enumerate_reduced_forms()
                .contains(&reduced_form)
        );
    }
}

#[test]
fn gp_split_ideal_form_fixtures_send_conjugates_to_inverse_classes() {
    let fixtures = gp_ideal_form_fixtures();
    for (left, right) in [
        (fixtures[0], fixtures[1]),
        (fixtures[2], fixtures[3]),
        (fixtures[4], fixtures[5]),
        (fixtures[6], fixtures[7]),
    ] {
        let left_ideal = left.ideal();
        let right_ideal = right.ideal();
        let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(left.discriminant))
            .expect("GP fixture discriminants should define class groups");

        assert_eq!(left.discriminant, right.discriminant);
        assert_eq!(left.ell, right.ell);
        assert_eq!(left_ideal.conjugate(), right_ideal);
        assert_eq!(right_ideal.conjugate(), left_ideal);

        let left_form = left.reduced_form();
        let right_form = right.reduced_form();
        assert_eq!(
            class_group
                .inverse(&left_form)
                .expect("GP fixture form should be a reduced class-group member"),
            right_form
        );
        assert_eq!(
            class_group
                .inverse(&right_form)
                .expect("GP fixture form should be a reduced class-group member"),
            left_form
        );
    }
}

#[test]
fn gp_ramified_ideal_form_fixture_is_fixed_by_conjugation() {
    let fixture = gp_ideal_form_fixtures()
        .into_iter()
        .find(|fixture| fixture.ramified)
        .expect("fixtures include one ramified case");
    let ideal = fixture.ideal();
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(fixture.discriminant))
        .expect("GP fixture discriminant should define a class group");
    let reduced_form = fixture.reduced_form();

    assert_eq!(ideal.conjugate(), ideal);
    assert_eq!(
        class_group
            .inverse(&reduced_form)
            .expect("GP fixture form should be a reduced class-group member"),
        reduced_form
    );
}
