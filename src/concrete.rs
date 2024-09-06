use crate::Ntt;
use concrete_ntt::prime64::{self, Plan};
use core::iter;
use rand_core::{OsRng, RngCore};

pub struct ConcreteNtt<const P: u64> {
    n: usize,
    plan: Plan,
}

impl<const P: u64> Ntt for ConcreteNtt<P> {
    type Elem = u64;

    fn new(n: usize) -> Self {
        Self {
            n,
            plan: prime64::Plan::try_new(n, P).unwrap(),
        }
    }

    fn rand(&self) -> Vec<Self::Elem> {
        iter::repeat_with(|| OsRng.next_u64() % P)
            .take(self.n)
            .collect()
    }

    fn forward<V: AsMut<[Self::Elem]>>(&self, mut a: V) -> V {
        self.plan.fwd(a.as_mut());
        a
    }

    fn backward<V: AsMut<[Self::Elem]>>(&self, mut a: V) -> V {
        self.plan.inv(a.as_mut());
        a
    }

    fn normalize<V: AsMut<[Self::Elem]>>(&self, mut a: V) -> V {
        self.plan.normalize(a.as_mut());
        a
    }
}
