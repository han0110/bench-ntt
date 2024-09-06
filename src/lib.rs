pub mod concrete;
pub mod plonky2;

pub trait Ntt {
    type Elem;

    fn new(n: usize) -> Self;

    fn rand(&self) -> Vec<Self::Elem>;

    fn forward<V: AsMut<[Self::Elem]>>(&self, a: V) -> V;

    fn backward<V: AsMut<[Self::Elem]>>(&self, a: V) -> V;

    fn normalize<V: AsMut<[Self::Elem]>>(&self, a: V) -> V;
}
