use ark_ed_on_bn254::Fq as F;

use super::copy_constraint::find_permutation;

pub fn gen_witness(x: F) -> (Vec<F>, Vec<F>, Vec<F>) {
    let a = vec![x, x * x, x * x * x, F::from(1), F::from(1), x * x * x + x];
    let b = vec![x, x, x, F::from(5), F::from(35), F::from(5)];
    let c = vec![
        x * x,
        x * x * x,
        x + x * x * x,
        F::from(5),
        F::from(35),
        F::from(35),
    ];

    (a, b, c)
}

pub fn is_satisfied_witness(a: Vec<F>, b: Vec<F>, c: Vec<F>) {
    assert!(a[0] * b[0] == c[0]);
    assert!(a[1] * b[1] == c[1]);
    assert!(a[2] + b[2] == c[2]);
    assert!(a[3] * b[3] == c[3]);
    assert!(a[4] * b[4] == c[4]);
    assert!(a[5] + b[5] == c[5]);
}

pub fn gen_copy_constraints() -> (Vec<F>, Vec<F>, Vec<F>, Vec<F>) {
    // copy constraints
    // a = [x , x*x, x*x*x, 1,  1, x*x*x + x]
    // b = [x , x, x, 5, 35, 5]
    // c = [x*x, x*x*x , x + x*x*x ,5, 35, 35]
    // inputs  = [x , x*x, x*x*x, 1,  1,
    //           x*x*x + x, x , x, x, 5, 35, 5
    //           x*x, x*x*x , x + x*x*x ,5, 35, 35]

    let copy_constraints = vec![
        F::from(8),
        F::from(12),
        F::from(13),
        F::from(3),
        F::from(4),
        F::from(14),
        F::from(0),
        F::from(6),
        F::from(7),
        F::from(15),
        F::from(17),
        F::from(9),
        F::from(1),
        F::from(2),
        F::from(5),
        F::from(11),
        F::from(10),
        F::from(16),
    ];

    let eval_domain: Vec<F> = (0..copy_constraints.len())
        .map(|i| F::from(i as i128))
        .collect();

    let x_a_prime = find_permutation(&copy_constraints[0..6], &eval_domain[0..6]);
    let x_b_prime = find_permutation(&copy_constraints[6..12], &eval_domain[6..12]);
    let x_c_prime = find_permutation(&copy_constraints[12..18], &eval_domain[12..18]);

    (x_a_prime, x_b_prime, x_c_prime, copy_constraints)
}
