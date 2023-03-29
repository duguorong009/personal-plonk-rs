use std::ops::Neg;

use ark_bn254::{Fr as ScalarField, G1Projective as G1, G2Projective as G2};
use ark_ec::Group;
use ark_ed_on_bn254::Fq as F;
use ark_ff::Zero;

pub fn powers_of_tau(secret: u128, depth: usize) -> (Vec<G1>, Vec<G2>) {
    let mut g_1 = vec![G1::generator()];
    let g_2 = vec![G2::generator(), G2::generator() * ScalarField::from(secret)];

    let mut sec = ScalarField::from(secret);
    for _ in 0..depth {
        g_1.push(G1::generator() * sec);
        sec = sec * ScalarField::from(secret);
    }

    (g_1, g_2)
}

pub fn poly_commit_g1(poly: &[F], g: &[G1]) -> G1 {
    let mut res = G1::zero();

    for (i, x) in poly.iter().enumerate() {
        // deal with negative numbers
        if (*x) < F::from(0) {
            let new_x = *x * F::from(-1);
            res += (g[i] * new_x).neg();
        } else {
            res += g[i] * x;
        }
    }

    res
}

pub fn poly_commit_g2(poly: &[F], g: &[G2]) -> G2 {
    let mut res = G2::zero();

    for (i, x) in poly.iter().enumerate() {
        // deal with negative numbers
        if (*x) < F::from(0) {
            let new_x = *x * F::from(-1);
            res += (g[i] * new_x).neg();
        } else {
            res += g[i] * x;
        }
    }

    res
}
