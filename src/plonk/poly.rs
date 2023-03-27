use ark_ed_on_bn254::Fq as F;

pub fn polynomial_eval(coef: &[F], x: F) -> F {
    let mut res = vec![];
    let mut power = F::from(1);
    for c in coef {
        res.push(*c * power);
        power *= x;
    }
    res.iter().sum()
}

pub fn polynomial_division() {
    todo!()
}

pub fn gen_poly() {
    todo!()
}
