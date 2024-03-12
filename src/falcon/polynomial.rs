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

use num_complex::Complex64;

use super::finite_field_element::FiniteFieldElem;
#[derive(Debug, Clone)]
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

    This will produce incorrect values in the case that the number of coefficients
    in p is not even.
     */
    pub fn split(p: &Polynomial) -> (Polynomial, Polynomial) {
        let length: usize = p.coefficients.len();

        let mut f0_coeffs: Vec<FiniteFieldElem> = Vec::new();
        let mut f1_coeffs: Vec<FiniteFieldElem> = Vec::new();

        for i in 0..(length / 2) {
            f0_coeffs.push(p.coefficients[2 * i].clone());
            f1_coeffs.push(p.coefficients[2 * i + 1].clone());
        }

        let f0: Polynomial = Polynomial::new(f0_coeffs);
        let f1: Polynomial = Polynomial::new(f1_coeffs);

        (f0, f1)
    }
    /*
        Splits a fft polynomial into two polynomials, one containing the even indices, and the other containing the odd.

        Floating point version of Polynomial::split(), used for FFT.
    */
    pub fn split_fp(f: &Vec<Complex64>) -> (Vec<Complex64>, Vec<Complex64>) {
        let mut f0: Vec<Complex64> = Vec::new();
        let mut f1: Vec<Complex64> = Vec::new();

        for (i, &value) in f.iter().enumerate() {
            if i % 2 == 0 {
                f0.push(value);
            } else {
                f1.push(value);
            }
        }

        (f0, f1)
    }

    /*
       Merges two polynomials into a single polynomial by interleaving their coefficients.

       Used after NTT computations to recombine the polynomial. Called in reverse order to that of
       split(). merge() is also called recursively in this manner.
    */
    pub fn merge(f_vec: Vec<Polynomial>) -> Polynomial {
        let f0 = f_vec.get(0).unwrap();
        let f1 = f_vec.get(1).unwrap();
        let length: usize = f0.coefficients.len() * 2;
        let mut f: Polynomial = Polynomial::new(vec![FiniteFieldElem::new(0); length]);

        for i in 0..length / 2 {
            f.coefficients[2 * i] = f0.coefficients[i];
            f.coefficients[2 * i + 1] = f1.coefficients[i];
        }

        f
    }

    /*
       Merges two fft polynomials into a single polynomial by interleaving their coefficients.

       Floating point version of Polynomial::merge(), used for FFT.
    */
    pub fn merge_fp(f0: Vec<f64>, f1: Vec<f64>) -> Vec<f64> {
        let length = f0.len() * 2;
        let mut f_merged = vec![0.0; length];

        for i in 0..f0.len() {
            f_merged[2 * i] = f0[i];
            f_merged[2 * i + 1] = f1[i];
        }

        f_merged
    }
}
