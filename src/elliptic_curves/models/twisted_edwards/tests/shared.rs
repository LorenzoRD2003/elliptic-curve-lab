use crate::elliptic_curves::TwistedEdwardsCurve;
use crate::fields::{Fp, traits::Field};

pub type F2 = Fp<2>;
pub type F5 = Fp<5>;
pub type F7 = Fp<7>;
pub type F13 = Fp<13>;

pub fn f5_curve() -> TwistedEdwardsCurve<F5> {
    TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular")
}

pub fn f13_curve() -> TwistedEdwardsCurve<F13> {
    TwistedEdwardsCurve::<F13>::new(F13::from_i64(3), F13::from_i64(5))
        .expect("sample twisted-Edwards curve should be non-singular")
}

pub fn f13_denominator_curve() -> TwistedEdwardsCurve<F13> {
    TwistedEdwardsCurve::<F13>::new(F13::from_i64(2), F13::one())
        .expect("sample twisted-Edwards curve should be non-singular")
}
