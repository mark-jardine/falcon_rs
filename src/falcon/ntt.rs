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

use std::vec;

use super::finite_field_element::FiniteFieldElem;
use super::ntt_consts::{INV_MOD_Q, ROOTS_DICT_ZQ};
use super::polynomial::Polynomial;

/*
    I2 - inverse of 2 mod Q
*/
const I2: u16 = 6145;

/*
    SQRT_1 - square root of (-1) mod q. ( found at ROOTS_DICT_ZQ.get(2)[0] )
*/
const SQRT_1: i32 = 1479;

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
    let i2_ffe: FiniteFieldElem = FiniteFieldElem::new(I2 as i32);

    let mut f0_ntt: Polynomial = Polynomial::new(vec![FiniteFieldElem::new(0); length / 2]);
    let mut f1_ntt: Polynomial = Polynomial::new(vec![FiniteFieldElem::new(0); length / 2]);

    for i in 0..length / 2 {
        //f0_ntt
        f0_ntt.coefficients[i].add(&f_ntt.coefficients[2 * i]);
        f0_ntt.coefficients[i].add(&f_ntt.coefficients[2 * i + 1]);
        f0_ntt.coefficients[i].mult(&i2_ffe);

        //f1_ntt
        f1_ntt.coefficients[i].sub(&f_ntt.coefficients[2 * i]);
        f1_ntt.coefficients[i].add(&f_ntt.coefficients[2 * i + 1]);
        f1_ntt.coefficients[i].mult(&i2_ffe);
        f1_ntt.coefficients[i].mult(&FiniteFieldElem::new(INV_MOD_Q[w[2 * i] as usize]));
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
fn merge_ntt(f_tup: (Polynomial, Polynomial)) -> Polynomial {
    let (f0_ntt, f1_ntt) = f_tup;
    let length: usize = f0_ntt.coefficients.len() * 2 as usize;
    let w = ROOTS_DICT_ZQ.get(&length).unwrap();
    let mut f_ntt: Polynomial = Polynomial::new(vec![FiniteFieldElem::new(0); length]);

    for i in 0..length / 2 {
        // Value of w[2*i]
        let w_at_2i = FiniteFieldElem::new(w.get(2 * i as usize).unwrap().to_owned());

        // f_ntt[2 * i + 0] = (f0_ntt[i] + w[2 * i] * f1_ntt[i]) % q
        let mut ffe: FiniteFieldElem = FiniteFieldElem::new(f0_ntt.coefficients[i].value as i32);
        ffe.add(&w_at_2i);
        ffe.mult(&f1_ntt.coefficients[i]);
        f_ntt.coefficients[2 * i] = ffe.clone();

        // f_ntt[2 * i + 1] = (f0_ntt[i] - w[2 * i] * f1_ntt[i]) % q
        ffe = FiniteFieldElem::new(f0_ntt.coefficients[i].value as i32);
        ffe.sub(&w_at_2i);
        ffe.mult(&f1_ntt.coefficients[i]);
        f_ntt.coefficients[2 * i + 1] = ffe.clone();
    }

    f_ntt
}

/*
    Compute the ntt of a polynomial

    Params:
        f - a polynomial in coefficient representation

    Returns:
        f_ntt - a polynomial in ntt representation
*/
fn ntt(f: &Polynomial) -> Polynomial {
    let len: usize = f.coefficients.len();
    let mut f_ntt: Polynomial = Polynomial::new(vec![]);

    if len > 2 {
        let (f0, f1) = Polynomial::split(&f);
        let f0_ntt: Polynomial = ntt(&f0);
        let f1_ntt: Polynomial = ntt(&f1);
        f_ntt = merge_ntt((f0_ntt, f1_ntt));
    } else if len == 2 {
        let sqrt_1_ffe = FiniteFieldElem::new(SQRT_1);

        f_ntt = Polynomial::new(vec![FiniteFieldElem::new(0); len]);

        // f_ntt[0] = (f[0] + sqr1 * f[1]) % q
        let mut ffe: FiniteFieldElem = f.coefficients[0];
        ffe.add(&sqrt_1_ffe);
        ffe.mult(&f.coefficients[1]);
        f_ntt.coefficients[0] = ffe.clone();

        // f_ntt[1] = (f[0] - sqr1 * f[1]) % q
        let mut ffe: FiniteFieldElem = f.coefficients[0];
        ffe.sub(&sqrt_1_ffe);
        ffe.mult(&f.coefficients[1]);
        f_ntt.coefficients[1] = ffe.clone();
    }

    f_ntt
}
