use crate::plonk::poly::polynomial_eval;
use ark_ed_on_bn254::Fq as F;
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

pub fn copy_constraint_simple(
    eval_domain: &[F],
    x_coef: &[F],
    y_coef: &[F],
    v1: F,
    v2: F,
) -> (Vec<F>, Vec<F>, Vec<F>, Vec<F>) {
    let mut p_x = vec![F::from(1)];
    let mut y = vec![];
    let mut rlc = vec![];
    let mut x = vec![];

    for i in 0..eval_domain.len() {
        x.push(polynomial_eval(&x_coef, eval_domain[i]));
        y.push(polynomial_eval(&y_coef, x[i]));

        rlc.push(v1 + x[i] + v2 * y[i]);
        p_x.push(p_x[i] * (v1 + x[i] + v2 * y[i]));
    }

    (x, y, p_x, rlc)
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
                &(&DensePolynomial::from_coefficients_slice(&[-x[k], F::one()])
                    / &DensePolynomial::from_coefficients_slice(&[fac])),
            );
        }
        p = p + pt;
    }

    p.coeffs
}

#[test]
fn test_lagrange() {
    let x = vec![
        F::from(0),
        F::from(1),
        F::from(2),
        F::from(3),
        F::from(4),
        F::from(5),
    ];
    let w = vec![
        F::from(0),
        F::from(1),
        F::from(2),
        F::from(3),
        F::from(4),
        F::from(5),
    ];
    let poly = lagrange(&x, &w);
    println!("{:?}", poly);
}
