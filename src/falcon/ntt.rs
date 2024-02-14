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

use std::error::Error;
use std::vec;

use super::finite_field_element::FiniteFieldElem;
use super::ntt_consts::{INV_MOD_Q, ROOTS_DICT_ZQ};
use super::polynomial::Polynomial;

const I2: u16 = 6145;

fn split_ntt(f_ntt: Polynomial) -> (Polynomial, Polynomial) {
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

    (Polynomial::new(vec![]), Polynomial::new(vec![]))
}
