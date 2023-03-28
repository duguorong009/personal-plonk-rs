use ark_ed_on_bn254::Fq as F;
use ark_ff::{BigInteger256, Field};

pub fn polynomial_eval_prime(
    coef: &[F],
    x: F,
    p: F,
    step_size: u128,   /* default value: 1 */
    start_power: u128, /* default value: 0 */
) -> F {
    let step_size = {
        let ss = F::from(step_size);
        let res: BigInteger256 = ss.into();
        res
    };

    let start_power = {
        let st = F::from(start_power);
        let res: BigInteger256 = st.into();
        res
    };

    let mut res = vec![];
    let mut power = x.pow(start_power);

    for i in coef {
        res.push(*i * power);
        power = power * x.pow(step_size);
    }

    // Just return "sum" since we do the computation in F_p field
    res.iter().sum::<F>()
}

pub fn fft(p: F, domain: &[F], poly: &[F]) -> Vec<F> {
    if poly.len() == 1 {
        return poly.to_vec();
    }

    let p_even: Vec<F> = (0..poly.len()).step_by(2).map(|i| poly[i]).collect();
    let p_odd: Vec<F> = (1..poly.len()).step_by(2).map(|i| poly[i]).collect();
    let domain_positive: Vec<F> = (0..domain.len()).step_by(2).map(|i| domain[i]).collect();

    let L = fft(p, &domain_positive, &p_even);
    let R = fft(p, &domain_positive, &p_odd);

    let mut p_x: Vec<F> = vec![];
    let mut p_x_minus_1: Vec<F> = vec![];

    for (i, (x, y)) in L.iter().zip(R.iter()).enumerate() {
        let y_times_root = *y * domain[i];

        p_x.push(*x + y_times_root);
        p_x_minus_1.push(*x - y_times_root);
    }

    let mut res = vec![];
    res.extend_from_slice(&p_x);
    res.extend_from_slice(&p_x_minus_1);

    res
}

pub fn ifft(p: F, domain: &[F], evaluation: &[F]) -> Vec<F> {
    let mut vals = fft(p, domain, evaluation);
    vals.reverse();
    vals.rotate_right(1);
    vals.iter()
        .map(|x| *x * F::from(vals.len() as u128).inverse().unwrap())
        .collect()
}
