mod plonk;
mod tests;

use ark_ed_on_bn254::Fq as F;

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

    // fn gen_witness(x: F) -> (Vec<F>, Vec<F>, Vec<F>) {
    //     let a = vec![x, x * x, x * x * x, F::from(1), F::from(1), x * x * x + x];
    //     let b = vec![x, x, x, F::from(5), F::from(35), F::from(5)];
    //     let c = vec![
    //         x * x,
    //         x * x * x,
    //         x + x * x * x,
    //         F::from(5),
    //         F::from(35),
    //         F::from(35),
    //     ];

    //     (a, b, c)
    // }

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

    use plonk::sample_problem::{gen_witness, is_satisfied_witness};

    // let (a, b, c) = gen_witness(F::from(1));

    // // Uncomment the next line and run
    // // is_satisfied_witness(a, b, c);

    // Reader should investigate why this fails and why the next passes
    let (a, b, c) = gen_witness(F::from(3));
    is_satisfied_witness(a, b, c);
}
