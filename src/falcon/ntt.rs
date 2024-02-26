/*
    Implementation of the Falcon Digital Signature Scheme.
    Copyright (C) 2024 Mark Jardine

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    Portions of this software are adapted from work by Thomas Prest, 2018,
    licensed under the MIT License. See the LICENSE_MIT file in this distribution
    for the full license text.
*/

use super::finite_field_element::{FiniteFieldElem, Q};
use super::ntt_consts::{INV_MOD_Q, ROOTS_DICT_ZQ};
use super::polynomial::Polynomial;
use std::vec;

#[derive(Debug)]
enum PolynomialError {
    ZeroLength,
}

/*
    I2 - inverse of 2 mod Q
*/
const I2: u16 = 6145;

/*
    SQRT_1 - square root of (-1) mod q. ( found at ROOTS_DICT_ZQ.get(2)[0] )
*/
const SQRT_1: u32 = 1479;

/*
    Split a polynomial f_ntt in ntt representation into
    two polynomials in ntt representation.

    Params:
        f_ntt: a polynomial in ntt representation

    Returns:
        A tuple of polynomials in ntt representation
*/
fn split_ntt(f_ntt: &Polynomial) -> (Polynomial, Polynomial) {
    let length: usize = f_ntt.coefficients.len();
    let w = ROOTS_DICT_ZQ.get(&length).unwrap();
    let mut f0_ntt: Polynomial = Polynomial::new(vec![FiniteFieldElem::new(0); length / 2]);
    let mut f1_ntt: Polynomial = Polynomial::new(vec![FiniteFieldElem::new(0); length / 2]);

    for i in 0..length / 2 {
        f0_ntt.coefficients[i].value = (I2 as u32
            * (f_ntt.coefficients[2 * i].value + f_ntt.coefficients[2 * i + 1].value))
            % Q as u32;
        let root = w.get(2 * i as usize).unwrap().to_owned() as usize;

        f1_ntt.coefficients[i].value = ((I2 as i64
            * (f_ntt.coefficients[2 * i].value as i64
                - f_ntt.coefficients[2 * i + 1].value as i64))
            * INV_MOD_Q[root] as i64)
            .rem_euclid(Q as i64) as u32;
    }

    (f0_ntt, f1_ntt)
}

/*
    Merge two polynomials in ntt representation into
    a single polynomial in ntt representation.

    Params:
        (f0_ntt, f1_ntt): a tuple of polynomials in ntt representation

    Returns:
        A polynomial in ntt representation
*/
fn merge_ntt(f_tup: (Polynomial, Polynomial)) -> Result<Polynomial, PolynomialError> {
    let (f0_ntt, f1_ntt) = f_tup;
    let length: usize = f0_ntt.coefficients.len() * 2 as usize;

    if length == 0 {
        return Err(PolynomialError::ZeroLength);
    }

    let w = ROOTS_DICT_ZQ.get(&length).unwrap();
    let mut f_ntt: Polynomial = Polynomial::new(vec![FiniteFieldElem::new(0); length]);

    for i in 0..length / 2 {
        // Value of w[2*i]
        let w_at_2i_i64 = w.get(2 * i).unwrap().to_owned() as i64;

        f_ntt.coefficients[2 * i].value = (f0_ntt.coefficients[i].value as i64
            + w_at_2i_i64 * f1_ntt.coefficients[i].value as i64)
            .rem_euclid(Q as i64) as u32;

        f_ntt.coefficients[2 * i + 1].value = (f0_ntt.coefficients[i].value as i64
            - w_at_2i_i64 * f1_ntt.coefficients[i].value as i64)
            .rem_euclid(Q as i64) as u32;
    }

    Ok(f_ntt)
}

/*
    Compute the ntt of a polynomial

    Params:
        f - a polynomial in coefficient representation

    Returns:
        f_ntt - a polynomial in ntt representation
*/
fn ntt(f: &Polynomial) -> Result<Polynomial, PolynomialError> {
    let len: usize = f.coefficients.len();
    let mut f_ntt: Polynomial = Polynomial::new(vec![]);

    if len > 2 {
        println!(
            "Before base case NTT: f.coefficients[0]: {}, f.coefficients[1]: {}",
            f.coefficients[0].value, f.coefficients[1].value
        );
        let (f0, f1) = Polynomial::split(&f);
        let f0_ntt: Polynomial = ntt(&f0)?;
        let f1_ntt: Polynomial = ntt(&f1)?;

        match merge_ntt((f0_ntt, f1_ntt)) {
            Ok(f) => {
                f_ntt = f;
            }
            Err(_) => {
                eprint!("Polynomial error in ntt().\n");
                return Err(PolynomialError::ZeroLength);
            }
        }
        println!(
            "After base case NTT: f_ntt.coefficients[0]: {}, f_ntt.coefficients[1]: {}",
            f_ntt.coefficients[0].value, f_ntt.coefficients[1].value
        );
    } else if len == 2 {
        f_ntt = Polynomial::new(vec![FiniteFieldElem::new(0); len]);

        f_ntt.coefficients[0].value = ((f.coefficients[0].value as i64
            + SQRT_1 as i64 * f.coefficients[1].value as i64)
            % Q as i64) as u32;

        f_ntt.coefficients[1].value = ((f.coefficients[0].value as i64
            - SQRT_1 as i64 * f.coefficients[1].value as i64)
            .rem_euclid(Q as i64)) as u32;
    }

    Ok(f_ntt)
}

fn inv_ntt(f_ntt: &Polynomial) -> Polynomial {
    let len: usize = f_ntt.coefficients.len();
    let mut f: Polynomial = Polynomial::new(vec![]);

    if len > 2 {
        println!(
            "Before base case iNTT: f_ntt.coefficients[0]: {}, f_ntt.coefficients[1]: {}",
            f_ntt.coefficients[0].value, f_ntt.coefficients[1].value
        );
        let (f0_ntt, f1_ntt) = split_ntt(&f_ntt);
        let f0 = inv_ntt(&f0_ntt);
        let f1 = inv_ntt(&f1_ntt);
        f = Polynomial::merge(vec![f0, f1]);
        println!(
            "After base case iNTT: f.coefficients[0]: {}, f.coefficients[1]: {}",
            f.coefficients[0].value, f.coefficients[1].value
        );
    } else if len == 2 {
        f = Polynomial::new(vec![FiniteFieldElem::new(0); len]);
        let inv_mod_q_at_1459 =
            FiniteFieldElem::new(*INV_MOD_Q.get(SQRT_1 as usize).unwrap() as u32);

        f.coefficients[0].value = (I2 as i64
            * (f_ntt.coefficients[0].value as i64 + f_ntt.coefficients[1].value as i64)
            % Q as i64) as u32;

        f.coefficients[1].value = ((I2 as i64
            * inv_mod_q_at_1459.value as i64
            * (f_ntt.coefficients[0].value as i64 - f_ntt.coefficients[1].value as i64)
                .rem_euclid(Q as i64))
        .rem_euclid(Q as i64)) as u32;
    }

    f
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_and_inv_ntt() {
        let f = Polynomial::new(vec![
            FiniteFieldElem::new(1),
            FiniteFieldElem::new(2),
            FiniteFieldElem::new(3),
            FiniteFieldElem::new(4),
        ]);

        // Compute the NTT of f
        let f_ntt = ntt(&f).expect("ntt() failed in test_ntt_and_inv_ntt");
        // println!("f after ntt computation: {f_ntt:?}");

        // Compute the inverse NTT of f_ntt
        let f_recovered = inv_ntt(&f_ntt);
        // println!("f after inv_ntt computation: {f_recovered:?}");

        // Check if f_recovered is equal to f
        assert_eq!(
            f.coefficients, f_recovered.coefficients,
            "ntt and inv_ntt are not inverses"
        );
    }

    #[test]
    fn test_ntt_zero_polynomial() {
        let zero_poly = Polynomial::new(vec![FiniteFieldElem::new(0); 4]);
        let ntt_poly = ntt(&zero_poly).expect("ntt() failed in test_ntt_zero_polynomial");

        assert!(
            ntt_poly.coefficients.iter().all(|&x| x.value == 0),
            "NTT of a zero polynomial should also be a zero polynomial."
        );
    }

    #[test]
    fn test_ntt_linearity() {
        let poly_a = Polynomial::new(vec![FiniteFieldElem::new(1), FiniteFieldElem::new(2)]);
        let poly_b = Polynomial::new(vec![FiniteFieldElem::new(3), FiniteFieldElem::new(4)]);

        // Compute sum of polynomials before NTT
        let mut sum_before_ntt = poly_a.clone();
        for (i, coeff) in poly_b.coefficients.iter().enumerate() {
            sum_before_ntt.coefficients[i].add(coeff);
        }

        // Compute NTT of the sum
        let ntt_sum = ntt(&sum_before_ntt).expect("ntt() failed in test_ntt_linearity");

        // Compute sum of NTTs
        let ntt_a = ntt(&poly_a).expect("ntt() failed in test_ntt_linearity");
        let ntt_b = ntt(&poly_b).expect("ntt() failed in test_ntt_linearity");
        let mut sum_of_ntts = ntt_a;
        for (i, coeff) in ntt_b.coefficients.iter().enumerate() {
            sum_of_ntts.coefficients[i].add(coeff);
        }

        assert_eq!(
            ntt_sum.coefficients, sum_of_ntts.coefficients,
            "NTT should be linear."
        );
    }

    #[test]
    fn test_ntt_random_polynomials() {
        use crate::falcon::finite_field_element::Q;
        use rand::{distributions::Uniform, Rng};

        let mut rng = rand::thread_rng();
        let dist = Uniform::from(0..Q as i64);

        for _ in 0..10 {
            let random_poly = Polynomial::new(
                (0..4)
                    .map(|_| FiniteFieldElem::new(rng.sample(&dist) as u32))
                    .collect(),
            );

            let ntt_poly = ntt(&random_poly).expect("ntt() failed in test_ntt_random_polynomials");
            let recovered_poly = inv_ntt(&ntt_poly);

            assert_eq!(
                random_poly.coefficients, recovered_poly.coefficients,
                "Random polynomial should be recovered accurately after NTT and inverse NTT."
            );
        }
    }

    #[test]
    fn test_degree_preservation() {
        let f = Polynomial::new(vec![
            FiniteFieldElem::new(1),
            FiniteFieldElem::new(3),
            FiniteFieldElem::new(5),
            FiniteFieldElem::new(7),
        ]);
        let f_ntt = ntt(&f).expect("ntt() failed in test_degree_preservation");

        let f_recovered = inv_ntt(&f_ntt);
        assert_eq!(
            f.coefficients.len(),
            f_recovered.coefficients.len(),
            "The degree of the polynomial should be preserved."
        );
    }
}
