use crate::plonk::poly::polynomial_eval;
use ark_ed_on_bn254::Fq as F;
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

pub fn copy_constraint_simple() {
    todo!()
}

pub fn find_permutation(copies: &[F], eval_domain: &[F]) -> Vec<F> {
    let perm = lagrange(eval_domain, copies);
    perm
}

fn lagrange<F>(x: &[F], w: &[F]) -> Vec<F>
where
    F: Field,
{
    let m = x.len();
    let mut p = DensePolynomial::from_coefficients_slice(&[F::zero()]);

    for j in 0..m {
        let mut pt = DensePolynomial::from_coefficients_slice(&[w[j]]);
        for k in 0..m {
            if k == j {
                continue;
            }
            let fac = x[j] - x[k];
            pt = pt.naive_mul(
                &(&DensePolynomial::from_coefficients_slice(&[F::one(), -x[k]])
                    / &DensePolynomial::from_coefficients_slice(&[fac])),
            );
        }
        p = p + pt;
    }
    p.coeffs.reverse();

    p.coeffs
}

#[test]
fn test_lagrange() {
    let x = vec![F::from(0), F::from(1), F::from(2)];
    let w = vec![F::from(4), F::from(5), F::from(6)];
    let poly = lagrange(&x, &w);
    println!("{:?}", poly);
}
