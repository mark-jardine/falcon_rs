use std::error::Error;

use num_bigint::{BigInt, Sign};

use super::{finite_field_element::Q, polynomial::Polynomial};

#[derive(Debug)]
pub enum NtruError {
    DivisionError,
    XgcdError,
    ReductionError,
    Other(String), // For miscellaneous errors
}

impl std::fmt::Display for NtruError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for NtruError {}

/*
   Compute the bitsize of an element of Z (not counting the sign).
   The bitsize is rounded to the next multiple of 8.
   This makes the function slightly imprecise, but faster to compute.
*/
pub fn bitsize(a: BigInt) -> u64 {
    // Will be a 61 1's, and the final three bits as 0's
    let bitmask = u64::MAX ^ 0b0111;
    // Add 7 to ensure rounding up to next multiple of 8
    // Then round down to a multiple of 8
    (a.bits() + 0b0111) & bitmask
    // E.g.
    // if a = 10
    // 10 + 7 & bitmask
    // 0000_1010 + 0000_0111 & (1111_1111..._1000)
    // = 0001_0001 [=17] & ..._1111_1000
    // = 0001_0000 [=16]
}

/*
    Compute the extended GCD of two integers b and n.
    Return d, u, v such that d = u * b + v * n, and d is the GCD of b, n.

    Params:
        - b: i64
        - n: i64

    Returns:
        - a tuple of three i64 values (d, u, v) such that d = u * b + v * n
*/
pub fn xgcd(mut b: BigInt, mut n: BigInt) -> (BigInt, BigInt, BigInt) {
    let (mut x0, mut x1, mut y0, mut y1) = (
        BigInt::from(1),
        BigInt::from(0),
        BigInt::from(0),
        BigInt::from(1),
    );

    while n != BigInt::from(0) {
        let q = b.clone() / n.clone();
        let n_tmp = n.clone();
        n = b % n;
        b = n_tmp;

        let x1_tmp = x1.clone();
        x1 = x0 - q.clone() * x1;
        x0 = x1_tmp;

        let y1_tmp = y1.clone();
        y1 = y0 - q * y1;
        y0 = y1_tmp;
    }

    (b, x0, y0)
}

pub fn ntru_solve(
    f: Polynomial<BigInt>,
    g: Polynomial<BigInt>,
) -> Result<(Polynomial<BigInt>, Polynomial<BigInt>), NtruError> {
    let len = f.coefficients.len();
    if len == 1 {
        let f0 = f.coefficients[0].clone();
        let g0 = g.coefficients[0].clone();
        let (d, u, v) = xgcd(f0, g0);
        if d != BigInt::from(1) {
            return Err(NtruError::XgcdError);
        } else {
            return Ok((
                Polynomial::new(vec![(v * -BigInt::from(Q))]),
                Polynomial::new(vec![BigInt::from(Q) * u]),
            ));
        }
    }

    Ok((f, g))
}
