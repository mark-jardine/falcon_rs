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
*/

/*
    A polynomial in coefficient format.
*/

use super::finite_field_element::FiniteFieldElem;
#[derive(Debug)]
pub struct Polynomial {
    pub coefficients: Vec<FiniteFieldElem>,
}

impl Polynomial {
    pub fn new(_coefficients: Vec<FiniteFieldElem>) -> Self {
        Polynomial {
            coefficients: _coefficients,
        }
    }
    /*
    Splits a polynomial into two polynomials, one containing the even indices, and the other containing the odd.

    Used in the initial step of NTT computation to increase efficiency by allowing recursive calls to smaller
    degree polynomials.

    todo: may require padding of zeros on p to ensure that the number of coefficients is even.
     */
    pub fn split(p: &Polynomial) -> (Polynomial, Polynomial) {
        let length: usize = p.coefficients.len();

        let mut f0_coeffs: Vec<FiniteFieldElem> = Vec::new();
        let mut f1_coeffs: Vec<FiniteFieldElem> = Vec::new();

        for i in 0..(length / 2) {
            f0_coeffs.push(p.coefficients[2 * i]);
            f1_coeffs.push(p.coefficients[2 * i + 1]);
        }

        let f0: Polynomial = Polynomial::new(f0_coeffs);
        let f1: Polynomial = Polynomial::new(f1_coeffs);

        (f0, f1)
    }

    /*
    Merges two polynomials into a single polynomial by interleaving their coefficients.

    Used after NTT computations to recombine the polynomial. Called in reverse order to that of
    split(). merge() is also called recursively in this manner.
     */
    pub fn merge(v: Vec<Polynomial>) -> Polynomial {
        todo!()
    }
}
