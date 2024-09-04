use core::any::type_name;
use criterion::{criterion_group, criterion_main, Criterion};
use ntt_bench::{
    concrete::ConcreteNtt,
    plonky2::{Plonky2Goldilocks, Plonky2Ntt},
    Ntt,
};

fn run<N: Ntt>(c: &mut Criterion, n: usize) {
    let ntt = N::new(n);
    let mut data = ntt.rand();
    c.bench_function(&format!("{}/{n}", type_name::<N>()), |b| {
        b.iter(|| ntt.forward(&mut data));
    });
}

fn goldilocks(c: &mut Criterion) {
    const P: u64 = 18446744069414584321;
    for n in [1024, 2048] {
        run::<Plonky2Ntt<Plonky2Goldilocks>>(c, n);
        run::<ConcreteNtt<P>>(c, n);
    }
}

criterion_group!(benches, goldilocks);
criterion_main!(benches);
