use criterion::{
    criterion_group, criterion_main, measurement::Measurement, BenchmarkGroup, Criterion,
};
use ntt_bench::{
    concrete::ConcreteNtt,
    plonky2::{Plonky2Goldilocks, Plonky2Ntt},
    Ntt,
};

fn drop<T>(_: T) {}

fn run_forward<N: Ntt>(c: &mut BenchmarkGroup<impl Measurement>, name: &str, n: usize) {
    let id = format!("{name}/{n}");
    let ntt = N::new(n);
    let mut a = ntt.rand();
    c.bench_function(id, |b| b.iter(|| drop(ntt.forward(&mut a))));
}

fn run_backward<N: Ntt>(c: &mut BenchmarkGroup<impl Measurement>, name: &str, n: usize) {
    let id = format!("{name}/{n}");
    let ntt = N::new(n);
    let mut a = ntt.rand();
    c.bench_function(id, |b| b.iter(|| drop(ntt.backward(&mut a))));
}

fn goldilocks(c: &mut Criterion) {
    const P: u64 = 18446744069414584321;
    {
        let g = &mut c.benchmark_group("forward");
        for n in [1024, 2048] {
            run_forward::<Plonky2Ntt<Plonky2Goldilocks>>(g, "plonky2", n);
            run_forward::<ConcreteNtt<P>>(g, "concrete", n);
        }
    }
    {
        let g = &mut c.benchmark_group("backward");
        for n in [1024, 2048] {
            run_backward::<Plonky2Ntt<Plonky2Goldilocks>>(g, "plonky2", n);
            run_backward::<ConcreteNtt<P>>(g, "concrete", n);
        }
    }
}

criterion_group!(benches, goldilocks);
criterion_main!(benches);
