mod plonk;
mod tests;

use ark_ed_on_bn254::Fq as F;
use itertools::izip;

fn main() {
    println!("Hello, world!");

    // What does our ZKP look like
    //
    // We have a bunch of stuff
    //  1. A witness which is information that lets us make the proof in our case this is knowledge of a variable x
    //      that satisfies our equation. It also includes the intermediate values that are zkp uses. This is secret to
    //      the user they want to prove this information.
    //  2. A set of gate constraints. Basically, that all the multiplications and additions we do are correct.
    //  3. A copy constraint check
    //  4. Input/output checks

    // Plonk Tutorial
    // This tutorial is based upon https://www.vitalik.ca/general/2019/09/22/plonk.html it expands upon the ideas
    // described there and teaches the user to build their own plonk implementation in rust.

    // This tutorial takes to forum of a series of challenges where the user eventually build a plonk implementation of a
    // single proof.

    // Plonk allows us to make arbitrary zero knowledge proofs. For the purposes of this tutorial we will prove that
    // we know an x such that P(x) = x^3 + x + 5 = 35  this is toy problem

    // [Diagram]

    // You can see we have two kinds of constraints - gate constraints and copy constraints. A constraint is like an "assert!"" from rust.
    // The program can only continue running if this assertion is true.

    // We will first handle the gate constraints and then tackle the copy constraints

    // Gen Witness
    // So lets first find a satisfying solution to the problem we are trying to make proofs about x^3 + x + 5 = 35

    // The variables a, b and c will be the checking of additions/multiplications operations. Where we define `a + b == c` or
    // `a * b == c`

    fn gen_witness(x: F) -> (Vec<F>, Vec<F>, Vec<F>) {
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

    // Right so now that we have a, b, c we are ready to test. `is_satisfied` tests that our witness matches
    // the constraints we are planning to add.

    // Basically a[0] * b[0] = c[0] check this is a multiplication
    // a[2] + b[2] = c[2] check this is addition

    // fn is_satisfied_witness(a: Vec<F>, b: Vec<F>, c: Vec<F>) {
    //     assert!(a[0] * b[0] == c[0]);
    //     assert!(a[1] * b[1] == c[1]);
    //     assert!(a[2] + b[2] == c[2]);
    //     assert!(a[3] * b[3] == c[3]);
    //     assert!(a[4] * b[4] == c[4]);
    //     assert!(a[5] + b[5] == c[5]);
    // }

    // use plonk::sample_problem::{gen_witness, is_satisfied_witness};

    // let (a, b, c) = gen_witness(F::from(1));

    // // Uncomment the next line and run
    // // is_satisfied_witness(a, b, c);

    // Reader should investigate why this fails and why the next passes
    let (a, b, c) = gen_witness(F::from(3));
    // is_satisfied_witness(a, b, c);

    // In plonk, everything is a polynomial

    // In the previous section we generated our witness. A witness is a valid solution to our constraints.
    // Where here our constraints are x^3 + x + 5 = 35

    // Next we want to define the actual constraints. They will be defined as a polynomial. Lets start out by
    // creating a function eval_poly which takes a polynomial and evaluates it at a given point. Take
    // the polynomial 1 + x + x^2 = y which is defined by this list [1, 1, 1] from lowest degree (ie starting at the 1 * x^0)
    // to (1 * x^2)

    // fn eval_poly(coef: &[F], x: F) -> F {
    //     let mut res = vec![];
    //     let mut power = F::from(1);
    //     for c in coef {
    //         res.push(*c * power);
    //         power *= x;
    //     }
    //     res.iter().sum()
    // }

    // assert!(eval_poly(&vec![F::from(1), F::from(1), F::from(1)], F::from(2)) == F::from(7));
    // assert!(
    //     eval_poly(
    //         &vec![F::from(-2), F::from(7), F::from(-5), F::from(1)],
    //         F::from(0)
    //     ) == F::from(-2)
    // );
    // assert!(
    //     eval_poly(
    //         &vec![F::from(-2), F::from(7), F::from(-5), F::from(1)],
    //         F::from(1)
    //     ) == F::from(1)
    // );
    // assert!(
    //     eval_poly(
    //         &vec![F::from(-2), F::from(7), F::from(-5), F::from(1)],
    //         F::from(2)
    //     ) == F::from(0)
    // );
    // assert!(
    //     eval_poly(
    //         &vec![F::from(-2), F::from(7), F::from(-5), F::from(1)],
    //         F::from(3)
    //     ) == F::from(1)
    // );

    // Okay now it seems out polynomial evaluations are working :)
    //
    // Our mul/add constraint is defined by (Q_l_i) * a_i + (Q_r_i) * b_i + (Q_o_i) * c_i + (Q_m_i) * a_i * b_i + Q_c_i = 0
    // we can use this to check additions and multiplications. Define the constant polynomial.

    // fn constaint_polynomial(
    //     Q_l_i: F,
    //     Q_r_i: F,
    //     Q_m_i: F,
    //     Q_o_i: F,
    //     Q_c_i: F,
    //     a_i: F,
    //     b_i: F,
    //     c_i: F,
    // ) -> bool {
    //     Q_l_i * a_i + Q_r_i * b_i + Q_o_i * c_i + Q_m_i * a_i * b_i + Q_c_i == F::from(0)
    // }

    // fn validate_native(
    //     Q_l: &[F],
    //     Q_r: &[F],
    //     Q_m: &[F],
    //     Q_o: &[F],
    //     Q_c: &[F],
    //     a: &[F],
    //     b: &[F],
    //     c: &[F],
    // ) -> bool {
    //     for (Q_l_i, Q_r_i, Q_m_i, Q_o_i, Q_c_i, a_i, b_i, c_i) in
    //         izip!(Q_l, Q_r, Q_m, Q_o, Q_c, a, b, c)
    //     {
    //         if constaint_polynomial(*Q_l_i, *Q_r_i, *Q_m_i, *Q_o_i, *Q_c_i, *a_i, *b_i, *c_i)
    //             == false
    //         {
    //             return false;
    //         }
    //     }
    //     true
    // }

    // fn test_addition() {
    //     // constraints
    //     let Q_l = vec![F::from(1)];
    //     let Q_r = vec![F::from(1)];
    //     let Q_m = vec![F::from(0)];
    //     let Q_o = vec![F::from(-1)];
    //     let Q_c = vec![F::from(0)];

    //     // witness
    //     let a = vec![F::from(0)];
    //     let b = vec![F::from(1)];
    //     let c = vec![F::from(1)];

    //     assert!(validate_native(&Q_l, &Q_r, &Q_m, &Q_o, &Q_c, &a, &b, &c));
    // }

    // fn test_mul() {
    //     // constraints
    //     let Q_l = vec![F::from(0)];
    //     let Q_r = vec![F::from(0)];
    //     let Q_m = vec![F::from(1)];
    //     let Q_o = vec![F::from(-1)];
    //     let Q_c = vec![F::from(0)];

    //     // witness
    //     let a = vec![F::from(1)];
    //     let b = vec![F::from(1)];
    //     let c = vec![F::from(1)];

    //     assert!(validate_native(&Q_l, &Q_r, &Q_m, &Q_o, &Q_c, &a, &b, &c));
    // }

    // fn test_constant() {
    //     // constraints
    //     let Q_l = vec![F::from(1)];
    //     let Q_r = vec![F::from(1)];
    //     let Q_m = vec![F::from(0)];
    //     let Q_o = vec![F::from(0)];
    //     let Q_c = vec![F::from(-10)];

    //     // witness
    //     let a = vec![F::from(10)];
    //     let b = vec![F::from(0)];
    //     let c = vec![F::from(10)];

    //     assert!(validate_native(&Q_l, &Q_r, &Q_m, &Q_o, &Q_c, &a, &b, &c));
    // }

    // test_addition();
    // test_mul();
    // test_constant();

    // Okay so now we are doing multiplications and additions. We can validate manually  that all of these
    // are being done correctly. Right so this is working we can make constraints. So lets make all the constraints
    // for our system. First lets make some helpers that drop the constraints where we need then make sure they
    // pass the tests.

    // use plonk::constraint::{add_add_constraint, add_constant_constraint, add_mul_constraint};

    // Okay now lets add the actual constraints. By setting Q_l, Q_r, Q_m, Q_o and Q_c such that it evaluates to a multiplication
    // constraint at Q[0] and an addition at Q[2].

    // fn gen_constraints() -> (Vec<F>, Vec<F>, Vec<F>, Vec<F>, Vec<F>) {
    //     // Prove that I know an X such that X * x * x + x + 5 == 35

    //     // Init constraints
    //     let mut Q_l = vec![];
    //     let mut Q_r = vec![];
    //     let mut Q_m = vec![];
    //     let mut Q_o = vec![];
    //     let mut Q_c = vec![];

    //     // set constraints
    //     add_mul_constraint(&mut Q_l, &mut Q_r, &mut Q_m, &mut Q_o, &mut Q_c);
    //     add_mul_constraint(&mut Q_l, &mut Q_r, &mut Q_m, &mut Q_o, &mut Q_c);
    //     add_add_constraint(&mut Q_l, &mut Q_r, &mut Q_m, &mut Q_o, &mut Q_c);
    //     add_constant_constraint(&mut Q_l, &mut Q_r, &mut Q_m, &mut Q_o, &mut Q_c, F::from(5));
    //     add_constant_constraint(
    //         &mut Q_l,
    //         &mut Q_r,
    //         &mut Q_m,
    //         &mut Q_o,
    //         &mut Q_c,
    //         F::from(35),
    //     );
    //     add_add_constraint(&mut Q_l, &mut Q_r, &mut Q_m, &mut Q_o, &mut Q_c);

    //     (Q_l, Q_r, Q_m, Q_o, Q_c)
    // }

    // Copy constraints

    // So at the moment the system is not secure. Basically, we are checking that the variables at location
    //  1. `a[0] * b[0] == c[0]`
    //  2. `a[1] * b[1] == c[1]`
    //  3. `a[2] + b[2] == c[2]`

    // But we are just hoping that `a[1] == c[0]` we need to add constraints to make sure that we copy c[0] to a[1] these are called
    // copy constraints you may also have heard of them referred to as permutation arguments.

    // Th naive thing to do is to do these check manually. Basically, make sure that each variable is equal to the other. The problem with this is
    // that it means that we need to check every variable which breaks privacy and succintness. Instead, we will find a way to do this check using
    // polynomials.

    // Right now we have our witness which is `witness = a + b + c` and we want to prove that the value at
    // `witness[0] == witness[8] == witness[7] == witness[6]` all of these values corresponding to our initial x.

    // So now we need to make 3 polynomials the first to return the index of witness we want to look up, `witness_x_1`.
    // The second witness `witness_x_2` to return the index after the permutation has been applied.
    // And the third `witness_y` which returns the actual value of that witness at a given index.

    // Write code that returns all of these.
    use crate::plonk::copy_constraint::find_permutation;
    use crate::plonk::poly::polynomial_eval;

    let mut witness = vec![];
    witness.extend_from_slice(&a);
    witness.extend_from_slice(&b);
    witness.extend_from_slice(&c);

    let eval_domain: Vec<F> = (0..witness.len()).map(|i| F::from(i as i128)).collect();
    let witness_x_a = find_permutation(
        &(0..a.len()).map(|i| F::from(i as i128)).collect::<Vec<F>>(),
        &(0..a.len()).map(|i| F::from(i as i128)).collect::<Vec<F>>(),
    );
    let witness_x_b = find_permutation(
        &(a.len()..(2 * b.len()))
            .map(|i| F::from(i as i128))
            .collect::<Vec<F>>(),
        &(a.len()..(2 * b.len()))
            .map(|i| F::from(i as i128))
            .collect::<Vec<F>>(),
    );
    let witness_x_c = find_permutation(
        &((2 * b.len())..(3 * c.len()))
            .map(|i| F::from(i as i128))
            .collect::<Vec<F>>(),
        &((2 * a.len())..(3 * a.len()))
            .map(|i| F::from(i as i128))
            .collect::<Vec<F>>(),
    );

    let witness_y = find_permutation(&witness, &eval_domain);
    for (i, val) in witness.iter().enumerate() {
        assert!(*val == polynomial_eval(&witness_y, F::from(i as i128)));
    }
}
