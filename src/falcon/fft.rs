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

use super::fft_consts::ROOTS_DICT;
use super::polynomial::Polynomial;
use num_complex::Complex64;
use std::vec;

/*
    Split a polynomial f into two polynomials

    Params:
        - f: a polynomial in FFT representation

    Returns:
        - a Tuple, containing two polynomials in FFT representation

    Reference: Algorithm 1 splitfft(FFT(f)) in Falcon's specification
*/
fn split_fft(f: Vec<Complex64>) -> (Vec<Complex64>, Vec<Complex64>) {
    let length = f.len();

    // Precomputed roots of unity associated with FFT of size `length`
    let w = ROOTS_DICT.get(&length).unwrap();

    let mut f0 = vec![Complex64::new(0.0, 0.0); length / 2];
    let mut f1 = vec![Complex64::new(0.0, 0.0); length / 2];

    for i in 0..(length / 2) {
        f0[i] = 0.5 * (f[2 * i] + f[2 * i + 1]);
        f1[i] = 0.5 * (f[2 * i] - f[2 * i + 1]) * w[2 * i].conj();
    }

    (f0, f1)
}
/*
    Merge two or three polynomials into a single polynomial f.

    Params:
        f_vec_fft: a vector of polynomials in FFT representation

    Reference: algorithm 2 (mergefft_2) in Falcon's specification.
*/
fn merge_fft(f_vec_fft: Vec<Vec<Complex64>>) -> Vec<Complex64> {
    let f0_fft: Vec<Complex64> = f_vec_fft[0].clone();
    let f1_fft: Vec<Complex64> = f_vec_fft[1].clone();

    let length: usize = 2 * f0_fft.len();
    let w: &Vec<Complex64> = ROOTS_DICT.get(&length).unwrap();
    let mut f_fft: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); length];

    for i in 0..(length / 2) {
        f_fft[2 * i] = f0_fft[i] + w[2 * i] * f1_fft[i];
        f_fft[2 * i + 1] = f0_fft[i] - w[2 * i] * f1_fft[i]
    }

    f_fft
}

/*
    Compute the FFT of a polynomial mod (x^n + 1)

    Params:
        - f: A polynomial in coefficient representation

    Returns:
        - A polynomial in FFT representation
*/
fn fft(f: Polynomial) -> Vec<Complex64> {
    let length: usize = f.coefficients.len();
    let mut f_fft: Vec<Complex64> = vec![];

    if length > 2 {
        let (f0, f1) = Polynomial::split(&f);
        let f0_fft = fft(f0);
        let f1_fft = fft(f1);
        f_fft = merge_fft(vec![f0_fft, f1_fft]);
    } else if length == 2 {
        f_fft = vec![Complex64::new(0.0, 0.0); length];
        f_fft[0] = Complex64::new(f.coefficients[0] as f64, 1.0 * f.coefficients[1] as f64);
        f_fft[1] = Complex64::new(f.coefficients[0] as f64, -1.0 * f.coefficients[1] as f64);
    }

    f_fft
}
