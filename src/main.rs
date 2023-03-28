mod plonk;
mod tests;

use ark_ed_on_bn254::Fq as F;
use ark_ff::{BigInteger, BigInteger256, Field, PrimeField};
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

    // The test code below checks that witness_y returns the same results when we use the permutated indexes or
    // the non permuated version. This means that each value there matches.

    // Okay now lets rearrange it so that the values get swapped when they match
    use crate::plonk::sample_problem::gen_copy_constraints;

    let (witness_x_a_perm, witness_x_b_perm, witness_x_c_perm, copy_constraints) =
        gen_copy_constraints();

    for i in 0..a.len() {
        assert!(
            polynomial_eval(
                &witness_y,
                polynomial_eval(&witness_x_a, F::from(i as i128))
            ) == polynomial_eval(
                &witness_y,
                polynomial_eval(&witness_x_a_perm, F::from(i as i128))
            )
        );
    }

    for i in a.len()..(2 * a.len()) {
        assert!(
            polynomial_eval(
                &witness_y,
                polynomial_eval(&witness_x_b, F::from(i as i128))
            ) == polynomial_eval(
                &witness_y,
                polynomial_eval(&witness_x_b_perm, F::from(i as i128))
            )
        );
    }

    for i in (2 * a.len())..(3 * a.len()) {
        assert!(
            polynomial_eval(
                &witness_y,
                polynomial_eval(&witness_x_c, F::from(i as i128))
            ) == polynomial_eval(
                &witness_y,
                polynomial_eval(&witness_x_c_perm, F::from(i as i128))
            )
        );
    }

    // So now we have a way of checking permutations with polynomials. But we still need to check every variable which means we
    // have not really gained anything. So next we will embed these three polynomials in a third such that we can check batches of
    // permutations at once.

    // To do this we take a random linear combination of witness_x_1 and witness_y.

    //Then we calculate P(x) where P(0) = 1 and

    // Then we do the same for witness_x_2 calculating P_2(x). Because v1 and v2 are random numbers we know that
    // P_1(i) == P_2(i) if and only if witness_Y gives the same results when evaluated on witness_x_1(0:i) and witness_x_2(0:i)

    use crate::plonk::copy_constraint::copy_constraint_simple;

    // We have to generate v1 and v2 after a, b and c have been fixed.
    let v1 = F::from(6263831568402553244_i128); // hash(str(a + b + c))
    let v2 = F::from(8550125927969859821_i128); // hash(str(c + b + a))

    let eval_domain: Vec<F> = (0..3 * a.len()).map(|i| F::from(i as i128)).collect();

    let (x, y, px_a, rlc) = copy_constraint_simple(
        &(0..a.len()).map(|i| F::from(i as i128)).collect::<Vec<F>>(),
        &witness_x_a,
        &witness_y,
        v1,
        v2,
    );

    let (x, y, px_b, rlc) = copy_constraint_simple(
        &(a.len()..2 * a.len())
            .map(|i| F::from(i as i128))
            .collect::<Vec<F>>(),
        &witness_x_b,
        &witness_y,
        v1,
        v2,
    );

    let (x, y, px_c, rlc) = copy_constraint_simple(
        &(2 * a.len()..3 * a.len())
            .map(|i| F::from(i as i128))
            .collect::<Vec<F>>(),
        &witness_x_c,
        &witness_y,
        v1,
        v2,
    );

    // calculate permutated polynomial
    let (x_1, y_1, px_a_prime, rlc_1) = copy_constraint_simple(
        &(0..a.len()).map(|i| F::from(i as i128)).collect::<Vec<F>>(),
        &witness_x_a_perm,
        &witness_y,
        v1,
        v2,
    );

    let (x_1, y_1, px_b_prime, rlc_1) = copy_constraint_simple(
        &(a.len()..2 * a.len())
            .map(|i| F::from(i as i128))
            .collect::<Vec<F>>(),
        &witness_x_b_perm,
        &witness_y,
        v1,
        v2,
    );

    let (x_1, y_1, px_c_prime, rlc_1) = copy_constraint_simple(
        &(2 * a.len()..3 * a.len())
            .map(|i| F::from(i as i128))
            .collect::<Vec<F>>(),
        &witness_x_c_perm,
        &witness_y,
        v1,
        v2,
    );

    assert!(
        px_a[px_a.len() - 1] * px_b[px_b.len() - 1] * px_c[px_c.len() - 1]
            == px_a_prime[px_a_prime.len() - 1]
                * px_b_prime[px_b_prime.len() - 1]
                * px_c_prime[px_c_prime.len() - 1]
    );
    assert!(px_a[0] == px_a_prime[0]);
    assert!(px_b[0] == px_b_prime[0]);
    assert!(px_c[0] == px_c_prime[0]);
    assert!(px_a[0] == px_b[0]);
    assert!(px_b[0] == px_c[0]);
    assert!(px_a[0] == F::from(1));

    // So now we can evaluate many copy constraints by simply checking a single point. But the problem is that the verifier needs to
    // compute the Px_a ... Px_c_prime. We want to come up with a way so that they don't need to evaluate these instead letting the
    // prover produce an argument that they have evaluated them correctly and minimize the verifiers work. We will do that in after the
    // next section. In the next section we will make a quick fft sidetrack cos we need that to make a performant prover.

    // Part x: FFT

    // So you can see that it takes about 1.5 seconds for 100 points. In reality we will want to make proofs for systems that constain
    // orders of magnatudes more variables in less than that time. So lets use FFT to speed it up.

    // So this is based upon https://vitalik.ca/general/2019/05/12/fft.html which is good to read before you continue.

    // TODO: possibly break this into a seperate tutorial

    // Firstly fft stands for fast foruier transform. A foruier transform is basically evaluating a polynomial.
    // There are two ways to represent polynomials

    // 1. Is via coefficients [1, 0, 3] is 1 + x^2 where its represented by coefficients
    // 2. Is with evaluations calculate the fourier space(evaluation space) values.

    // use crate::plonk::poly::polynomial_eval;

    let poly = vec![F::from(1), F::from(0), F::from(3)];
    let mut fs = vec![];
    fs.push(polynomial_eval(&poly, F::from(0)));
    fs.push(polynomial_eval(&poly, F::from(1)));
    fs.push(polynomial_eval(&poly, F::from(2)));

    // println!("fs: {:?}", fs);

    // Both coefficient space and evaluation space uniquely identify the polynomial as long as it has been evaluated at a few
    // positions. So what we did by evaluating the polynomial is a fourier transform. It is also possible to do an inverse foruier
    // transform by basically interpolating the polynomial to find the coefficient form.
    use crate::plonk::copy_constraint::lagrange;

    let x = vec![F::from(0), F::from(1), F::from(2)];
    let res = lagrange(&x, &fs);
    assert!(res == poly);

    // Okay so we have gone to foruier space and back to coefficient space. But why?

    // Well turns out in evaluation space it is easier to do things like multiplicaion and division. So lets do that now. Lets take the polynomial
    // 1 + 3 * x^2 and multiply it by itself. Lets do it both ways in fourier space and in coordinate space.

    let mut res: Vec<F> = (0..(poly.len().pow(2))).map(|_| F::from(0)).collect();
    let expected_result = vec![
        F::from(1),
        F::from(0),
        F::from(6),
        F::from(0),
        F::from(9),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
    ];
    for (i, coef1) in poly.iter().enumerate() {
        for (j, coef2) in poly.iter().enumerate() {
            res[i + j] += coef1.clone() * coef2.clone();
        }
    }
    assert!(res == expected_result);

    // Okay so this took len(poly) ^ 2 operations. Lets try the fft version now.
    let poly = vec![F::from(1), F::from(0), F::from(3)];
    let mut fs = vec![];
    for i in 0..9 {
        fs.push(polynomial_eval(&poly, F::from(i as i128)));
    }

    let fs_res: Vec<F> = fs
        .iter()
        .zip(fs.iter())
        .map(|(x, y)| x.clone() * y.clone())
        .collect();

    let x: Vec<F> = (0..9).map(|i| F::from(i as i128)).collect();
    let res = lagrange(&x, &fs_res);
    // assert!(res == expected_result);  // since we use DensePolynomial of arkworks, it strips the zeros at the end

    // Okay so this tool len(poly)*2 operations to do the multiplicaions. But there is still a problem. Because the transform and inverse
    // both cost more than the saving in operations we need to find a faster way to do this fourier transform. That is where the fast in
    // fast fourier transform comes in.

    // Okay so firstly lets evaluate a polynomial over prime feild. You just have to take the result % p.

    use crate::plonk::fft::polynomial_eval_prime;
    // // this evaluations 3 + x**2 at position 0 % 5
    // assert!(polynomial_eval_prime(&vec![3, 0, 1], 0, 5, 1, 0) == F::from(3));
    // // this evaluations 3 + x**2 at position 1 % 5
    // assert!(polynomial_eval_prime(&vec![3, 0, 1], 1, 5, 1, 0) == F::from(4));

    // // this evaluations 3 + x**2 at position 2 % 5
    // assert!(polynomial_eval_prime(&vec![3, 0, 1], 2, 5, 1, 0) == F::from(2));

    // Every prime feild has something called roots of unity. Which is a number x such that x^n == 1.

    // This evaluation 3 + x^2 at position 1
    fn roots_of_unity(order: usize) -> Vec<F> {
        let a = F::from(5);
        let modulus = <F as PrimeField>::MODULUS;
        let p: F = F::from_le_bytes_mod_order(&modulus.to_bytes_le());

        let mut res = vec![];
        for i in 0..order {
            let exp: BigInteger256 =
                (F::from(i as i128) * (p - F::from(1)) / F::from(order as i128)).into();
            res.push(a.pow(exp));
        }
        res
    }
    let roots_2 = roots_of_unity(2);

    // So lets take these roots and square them see what happens.
    // assert!(roots_2[0].pow([2]) == roots_2[0]);
    assert!(roots_2[1].pow([2]) == F::from(1));
    assert!(roots_2[1].pow([3]) == roots_2[1]);
    assert!(roots_2[1].pow([4]) == F::from(1));

    // We see that they stay the same. Let's say that we have a polynomial above 3 + x + x^2 + x^3 + x^4. If we want to evaluate it
    // at a point x^n = 1 where n = 2.

    // Let's say that we want to evaluate this polynomial at two points the native thing to do is to hen we can simply the above equation
    // because we know that x, x^3 == roots_2[1] and

    // Now, lets use the roots of unity to accelerate the evaluation of the polynomial above.

    // So we can save doing the squaring above because x^2 == x^4 for these roots. That lets us do the evaluation by summing 3 +
    // roots[0] + roots[0] and 3 + roots[1] + roots[1].

    let modulus = <F as PrimeField>::MODULUS;
    let p: F = F::from_le_bytes_mod_order(&modulus.to_bytes_le());
    let eval1 = F::from(3) + roots_2[0] + roots_2[0] + roots_2[0] + roots_2[0];
    let eval2 = F::from(3) + roots_2[1] + roots_2[0] + roots_2[1] + roots_2[0];

    assert!(
        polynomial_eval_prime(
            &[F::from(3), F::from(1), F::from(1), F::from(1), F::from(1)],
            roots_2[0],
            p,
            1,
            0,
        ) == eval1
    );
    assert!(
        polynomial_eval_prime(
            &[F::from(3), F::from(1), F::from(1), F::from(1), F::from(1)],
            roots_2[1],
            p,
            1,
            0
        ) == eval2
    );
    println!("{}", eval2);

    // So the last part is a bit strange. Since eval2 == 3 even tho roots[1] = a pretty big number. The reason for this is that
    // roots[0] == -roots[1] read about negative numbers in finite fields https://vitalik.ca/general/2017/11/22/starks_part_2.html

    // So because a negative number squared is a positive number we only have to worry about half the domain for even powers. Check out
    // the rest of https://vitalik.ca/general/2019/05/12/fft.html to fill in the details.

    // Then write an fft that validates the testcase below
    use crate::plonk::fft::fft;

    let domain = roots_of_unity(8);
    let poly = vec![
        F::from(3),
        F::from(1),
        F::from(4),
        F::from(1),
        F::from(5),
        F::from(9),
        F::from(2),
        F::from(6),
    ];
    let p_x = fft(p, &domain, &poly);

    let mut result = vec![];
    for x in domain {
        result.push(polynomial_eval_prime(&poly, x, p, 1, 0));
    }
    assert!(p_x == result);

    // And the same for ifft
    use crate::plonk::fft::ifft;

    let domain = roots_of_unity(8);
    let poly = vec![
        F::from(3),
        F::from(1),
        F::from(4),
        F::from(1),
        F::from(5),
        F::from(9),
        F::from(2),
        F::from(6),
    ];
    let p_x = fft(p, &domain, &poly);

    let result = ifft(p, &domain, &p_x);
    assert!(result == poly);

    // Okay thats fft done. We can do fast evaluations of polynomials. Do multiplicaion and stuff in fourier space and then convert back
    // to coefficient space.

    // TODO: do this better

    // Part x: Division of polynomials

    // Okay, the next thing we need to do is find how to divide polynomials efficiently. We can work out a basic algorithm that does this but
    // lets do it in fourier space since that will work out as being a bit faster.

    // Take the polynomial 1 + x and square it using fft.

    let poly1 = vec![
        F::from(1),
        F::from(1),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
    ];
    let poly2 = vec![
        F::from(1),
        F::from(1),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
    ];

    let domain = roots_of_unity(8);

    let poly1_fs = fft(p, &domain, &poly1);
    let poly2_fs = fft(p, &domain, &poly2);

    let res: Vec<F> = poly1_fs
        .iter()
        .zip(poly2_fs.iter())
        .map(|(x, y)| *x * *y)
        .collect();
    let res = ifft(p, &domain, &res);
    assert!(
        res == vec![
            F::from(1),
            F::from(2),
            F::from(1),
            F::from(0),
            F::from(0),
            F::from(0),
            F::from(0),
            F::from(0)
        ]
    );

    // Good now divide hint: remember to use field.FQ for it the division

    let poly1 = vec![F::from(1), F::from(2), F::from(1), F::from(0)];
    let poly2 = vec![F::from(1), F::from(1), F::from(0), F::from(0)];

    let domain = roots_of_unity(4);

    let poly1_fs = fft(p, &domain, &poly1);
    let poly2_fs = fft(p, &domain, &poly2);
    let res: Vec<F> = poly1_fs
        .iter()
        .zip(poly2_fs.iter())
        .map(|(x, y)| if *y == F::from(0) { *x * *y } else { *x / *y })
        .collect();

    let res = ifft(p, &domain, &res);
    assert!(res == vec![F::from(1), F::from(1), F::from(0), F::from(0)]);

    // Okay we can mutiply and divide in fft form which is cool. But what happens when we try and divide a polynomial by one that does not
    // have a root?

    let poly1 = vec![F::from(1), F::from(1), F::from(1), F::from(0)];
    let poly2 = vec![F::from(1), F::from(1), F::from(0), F::from(0)];

    let domain = roots_of_unity(8);

    let poly1_fs = fft(p, &domain, &poly1);
    let poly2_fs = fft(p, &domain, &poly2);
    let res: Vec<F> = poly1_fs
        .iter()
        .zip(poly2_fs.iter())
        .map(|(x, y)| if *y == F::from(0) { *x * *y } else { *x / *y })
        .collect();

    let res = ifft(p, &domain, &res);
    println!("{:?}", res);

    // TODO: find out the meaning of this polynomial

    // Okay so now we are able to quickly divide polynomials using fft what do we want this for?

    // Earlier we had polynomials that the verifier needed to multiply together in order to make sure they matched the data the prover sent.
    // This work was square in the size of the polynomials which meant that our proof were not succinct. In order to make our proofs succinct
    // we need to turn the verification into a bunch of polynomial evaluations.

    // Part x: is zero check

    // So we are working with polynomials and we want to be sure that a polynomial equals zero everywhere inside a domain that we care
    // about. So what we do is make it so that the roots of an equation are equal to zero everywhere we care about.

    // So we have a polynomial f(x) = and we want to prove that it is zero at a bunch of places we care about say (x=1, x=2)

    //This is easy to do all we need to do is come up with a list of places that we care about (1,2,3) and rewrite our polynomial as the product of
    // z(x) * h(x) == f(x) where z(x) = (x - 1)(x - 2)

    // So we have the polynomial f(x) = x^4 - 10x^3 + 35x^2 - 50x + 24 and we have the polynomial
    // (x - 2)(x - 3) = x^2 - 5x + 6

    // Use the polynomial division from above to calculate h(x)

    let fx = vec![
        F::from(24),
        F::from(-50),
        F::from(35),
        F::from(-10),
        F::from(1),
        F::from(0),
        F::from(0),
        F::from(0),
    ];
    let zx = vec![
        F::from(6),
        F::from(-5),
        F::from(1),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
        F::from(0),
    ];

    let domain = roots_of_unity(8);

    let fx_fs = fft(p, &domain, &fx);
    let zx_fs = fft(p, &domain, &zx);

    let res: Vec<F> = fx_fs
        .iter()
        .zip(zx_fs.iter())
        .map(|(x, y)| if *y == F::from(0) { *x * *y } else { *x / *y })
        .collect();
    let res = ifft(p, &domain, &res);

    // Convert to negative representation
    let hx: Vec<F> = res
        .iter()
        .map(|i| if *i > (p / F::from(2)) { -(p - *i) } else { *i })
        .collect();

    assert!(
        hx == vec![
            F::from(4),
            F::from(-5),
            F::from(1),
            F::from(0),
            F::from(0),
            F::from(0),
            F::from(0),
            F::from(0)
        ]
    );
    // [4, -5, 1] == (x - 1)(x - 4)

    // todo: for some reason when i divide by x-1 it fails check this out.

    // Okay so now we have a way of checking that a poylnomial == 0 at every point. That we define. Next what we want to do is remove the
    // requirement for the verifier to multiply two polynomials together and check that the results are equal.

    // So we avoid this by instead of checking the polynomial at every point we check it at a single random point. So the verifier has
    // polynomial

    // f(x) and h(x) * z(x) so what we do is generate a random number (rand) and evaluate f(x) and h(x), z(x) at rand then we
    // assert that f(rand) == h(rand) * z(rand)

    let rand = F::from(6201848214169226457_u128); // rand = hash(str(fx))
    assert!(
        polynomial_eval_prime(&fx, rand, p, 1, 0)
            == polynomial_eval_prime(&zx, rand, p, 1, 0)
                * polynomial_eval_prime(&hx, rand, p, 1, 0)
    );

    // Turns out that this is the case for the majority of points on fx , zx and hx. Becuase we generate the rand point to evaluate fx we will
    // only know that after it is created. So an attacker would have to spend a very long time to try and generate a random number that
    // passes the test. So the the thinking is that this is enough to secure plonk and only check a single point.
}
