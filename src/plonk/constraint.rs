use ark_ed_on_bn254::Fq as F;
use itertools::izip;

pub fn constaint_polynomial(
    Q_l_i: F,
    Q_r_i: F,
    Q_m_i: F,
    Q_o_i: F,
    Q_c_i: F,
    a_i: F,
    b_i: F,
    c_i: F,
) -> bool {
    Q_l_i * a_i + Q_r_i * b_i + Q_o_i * c_i + Q_m_i * a_i * b_i + Q_c_i == F::from(0)
}

pub fn validate_native(
    Q_l: &[F],
    Q_r: &[F],
    Q_m: &[F],
    Q_o: &[F],
    Q_c: &[F],
    a: &[F],
    b: &[F],
    c: &[F],
) -> bool {
    for (Q_l_i, Q_r_i, Q_m_i, Q_o_i, Q_c_i, a_i, b_i, c_i) in
        izip!(Q_l, Q_r, Q_m, Q_o, Q_c, a, b, c)
    {
        if constaint_polynomial(*Q_l_i, *Q_r_i, *Q_m_i, *Q_o_i, *Q_c_i, *a_i, *b_i, *c_i) == false {
            return false;
        }
    }
    true
}

pub fn add_add_constraint(
    Q_l: &mut Vec<F>,
    Q_r: &mut Vec<F>,
    Q_m: &mut Vec<F>,
    Q_o: &mut Vec<F>,
    Q_c: &mut Vec<F>,
) {
    Q_l.push(F::from(1));
    Q_r.push(F::from(1));
    Q_m.push(F::from(0));
    Q_o.push(F::from(-1));
    Q_c.push(F::from(0));
}

pub fn add_mul_constraint(
    Q_l: &mut Vec<F>,
    Q_r: &mut Vec<F>,
    Q_m: &mut Vec<F>,
    Q_o: &mut Vec<F>,
    Q_c: &mut Vec<F>,
) {
    Q_l.push(F::from(0));
    Q_r.push(F::from(0));
    Q_m.push(F::from(1));
    Q_o.push(F::from(-1));
    Q_c.push(F::from(0));
}

pub fn add_constant_constraint(
    Q_l: &mut Vec<F>,
    Q_r: &mut Vec<F>,
    Q_m: &mut Vec<F>,
    Q_o: &mut Vec<F>,
    Q_c: &mut Vec<F>,
    constant: F,
) {
    Q_l.push(F::from(0));
    Q_r.push(F::from(0));
    Q_m.push(F::from(1));
    Q_o.push(F::from(0));
    Q_c.push(F::from(-constant));
}
