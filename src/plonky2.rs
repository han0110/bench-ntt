use crate::Ntt;
use itertools::izip;
use plonky2_field::{packable::Packable, packed::PackedField, types::Field};
use unroll::unroll_for_loops;

pub use plonky2_field::goldilocks_field::GoldilocksField as Plonky2Goldilocks;

pub struct Plonky2Ntt<F: Field> {
    n: usize,
    log_n: usize,
    log_n_packed: usize,
    twiddles: Vec<F>,
}

impl<F: Field> Ntt for Plonky2Ntt<F> {
    type Elem = F;

    fn new(n: usize) -> Self {
        assert!(n.is_power_of_two());
        let log_n = n.ilog2() as _;
        Self {
            n,
            log_n,
            log_n_packed: log_n.saturating_sub(<F as Packable>::Packing::WIDTH.ilog2() as _),
            twiddles: bit_reverse(
                F::primitive_root_of_unity(log_n + 1)
                    .powers()
                    .take(n)
                    .collect(),
            ),
        }
    }

    fn rand(&self) -> Vec<Self::Elem> {
        F::rand_vec(self.n)
    }

    #[unroll_for_loops]
    fn forward(&self, a: &mut [Self::Elem]) {
        for layer in 0..self.log_n_packed {
            let a = <F as Packable>::Packing::pack_slice_mut(a);
            let (m, size) = (1 << layer, 1 << (self.log_n_packed - 1 - layer));
            izip!(a.chunks_exact_mut(2 * size), &self.twiddles[m..]).for_each(|(a, t)| {
                let (a, b) = a.split_at_mut(size);
                izip!(a, b).for_each(|(a, b)| dit(a, b, t));
            });
        }

        for i in 0..4 {
            let layer = self.log_n.saturating_sub(4 - i);
            if (self.log_n_packed..self.log_n).contains(&layer) {
                let (m, size) = (1 << layer, 1 << (3 - i));
                izip!(a.chunks_exact_mut(2 * size), &self.twiddles[m..]).for_each(|(a, t)| {
                    let (a, b) = a.split_at_mut(size);
                    izip!(a, b).for_each(|(a, b)| dit(a, b, t));
                });
            }
        }
    }
}

fn bit_reverse<T, V: AsMut<[T]>>(mut a: V) -> V {
    let n = a.as_mut().len();
    if n > 2 {
        assert!(n.is_power_of_two());
        let log_n = n.ilog2();
        for i in 0..n {
            let j = i.reverse_bits() >> (usize::BITS - log_n);
            if i < j {
                a.as_mut().swap(i, j)
            }
        }
    }
    a
}

#[inline(always)]
fn dit<F: Field, P: PackedField<Scalar = F>>(a: &mut P, b: &mut P, t: &F) {
    let bt = *b * *t;
    let c = *a + bt;
    let d = *a - bt;
    *a = c;
    *b = d;
}

#[cfg(test)]
mod test {
    use crate::{
        concrete::ConcreteNtt,
        plonky2::{Plonky2Goldilocks, Plonky2Ntt},
        Ntt,
    };

    #[test]
    fn forward() {
        const P: u64 = 18446744069414584321;
        for n in (4..12).map(|log_n| 1 << log_n) {
            let plonky2 = Plonky2Ntt::<Plonky2Goldilocks>::new(n);
            let concrete = ConcreteNtt::<P>::new(n);
            let mut a = plonky2.rand();
            let mut b = a.iter().map(|data| data.0).collect::<Vec<_>>();
            plonky2.forward(&mut a);
            concrete.forward(&mut b);
            assert_eq!(a.iter().map(|data| data.0).collect::<Vec<_>>(), b);
        }
    }
}
