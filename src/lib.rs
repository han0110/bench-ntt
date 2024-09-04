pub mod concrete;
pub mod plonky2;

pub trait Ntt {
    type Elem;

    fn new(n: usize) -> Self;

    fn rand(&self) -> Vec<Self::Elem>;

    fn forward(&self, a: &mut [Self::Elem]);
}
