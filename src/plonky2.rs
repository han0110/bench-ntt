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
    twiddle_invs: Vec<F>,
    n_inv: F,
}

impl<F: Field> Ntt for Plonky2Ntt<F> {
    type Elem = F;

    fn new(n: usize) -> Self {
        assert!(n.is_power_of_two());
        assert!(<F as Packable>::Packing::WIDTH.is_power_of_two());
        let log_n = n.ilog2() as _;
        let omega = F::primitive_root_of_unity(log_n + 1);
        Self {
            n,
            log_n,
            log_n_packed: log_n.saturating_sub(<F as Packable>::Packing::WIDTH.ilog2() as _),
            twiddles: bit_reverse(omega.powers().take(n).collect()),
            twiddle_invs: bit_reverse(omega.inverse().powers().take(n).collect()),
            n_inv: F::from_canonical_usize(n).inverse(),
        }
    }

    fn rand(&self) -> Vec<Self::Elem> {
        F::rand_vec(self.n)
    }

    #[unroll_for_loops]
    fn forward<V: AsMut<[Self::Elem]>>(&self, mut a: V) -> V {
        for layer in 0..self.log_n_packed {
            let a = <F as Packable>::Packing::pack_slice_mut(a.as_mut());
            let (m, size) = (1 << layer, 1 << (self.log_n_packed - 1 - layer));
            izip!(a.chunks_exact_mut(2 * size), &self.twiddles[m..]).for_each(|(a, t)| {
                let (a, b) = a.split_at_mut(size);
                izip!(a, b).for_each(|(a, b)| dit(a, b, t));
            });
        }
        for i in 0..4 {
            let layer = self.log_n.wrapping_sub(4 - i);
            if (self.log_n_packed..self.log_n).contains(&layer) {
                let a = a.as_mut();
                let (m, size) = (1 << layer, 1 << (3 - i));
                izip!(a.chunks_exact_mut(2 * size), &self.twiddles[m..]).for_each(|(a, t)| {
                    let (a, b) = a.split_at_mut(size);
                    izip!(a, b).for_each(|(a, b)| dit(a, b, t));
                });
            }
        }
        a
    }

    #[unroll_for_loops]
    fn backward<V: AsMut<[Self::Elem]>>(&self, mut a: V) -> V {
        for i in 0..4 {
            let layer = self.log_n.wrapping_sub(i + 1);
            if (self.log_n_packed..self.log_n).contains(&layer) {
                let a = a.as_mut();
                let (m, size) = (1 << layer, 1 << i);
                izip!(a.chunks_exact_mut(2 * size), &self.twiddle_invs[m..]).for_each(|(a, t)| {
                    let (a, b) = a.split_at_mut(size);
                    izip!(a, b).for_each(|(a, b)| dif(a, b, t));
                });
            }
        }
        for layer in (0..self.log_n_packed).rev() {
            let a = <F as Packable>::Packing::pack_slice_mut(a.as_mut());
            let (m, size) = (1 << layer, 1 << (self.log_n_packed - 1 - layer));
            izip!(a.chunks_exact_mut(2 * size), &self.twiddle_invs[m..]).for_each(|(a, t)| {
                let (a, b) = a.split_at_mut(size);
                izip!(a, b).for_each(|(a, b)| dif(a, b, t));
            });
        }
        a
    }

    fn normalize<V: AsMut<[Self::Elem]>>(&self, mut a: V) -> V {
        if self.n >= <F as Packable>::Packing::WIDTH {
            let a = <F as Packable>::Packing::pack_slice_mut(a.as_mut());
            a.iter_mut().for_each(|a| *a *= self.n_inv);
        } else {
            let a = a.as_mut();
            a.iter_mut().for_each(|a| *a *= self.n_inv);
        }
        a
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
    (*a, *b) = (*a + bt, *a - bt)
}

#[inline(always)]
fn dif<F: Field, P: PackedField<Scalar = F>>(a: &mut P, b: &mut P, t: &F) {
    (*b, *a) = ((*a - *b) * *t, *a + *b);
}

#[cfg(test)]
mod test {
    use crate::{
        plonky2::{bit_reverse, Plonky2Goldilocks, Plonky2Ntt},
        Ntt,
    };
    use plonky2_field::polynomial::PolynomialCoeffs;

    #[test]
    fn forward() {
        for n in (0..12).map(|log_n| 1 << log_n) {
            let ntt = Plonky2Ntt::<Plonky2Goldilocks>::new(n);
            let a = ntt.rand();
            let omega = ntt.twiddles[n / 2];
            assert_eq!(
                ntt.forward(a.clone()),
                bit_reverse(PolynomialCoeffs::new(a).coset_fft(omega).values)
            );
        }
    }

    #[test]
    fn round_trip() {
        for n in (0..12).map(|log_n| 1 << log_n) {
            let ntt = Plonky2Ntt::<Plonky2Goldilocks>::new(n);
            let a = ntt.rand();
            assert_eq!(ntt.normalize(ntt.backward(ntt.forward(a.clone()))), a);
        }
    }
}
