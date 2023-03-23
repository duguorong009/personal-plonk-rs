use ark_ed_on_bn254::Fq as F;

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
