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
use num_complex::Complex64;

/*
    Split a polynomial f into two polynomials

    Params:
        - f: a polynomial in FFT representation

    Returns:
        - a Tuple, containing two polynomials

    Reference: Algorithm 1 splitfft(FFT(f)) in Falcon specification
*/
fn split_fft(f: Vec<Complex64>) -> (Vec<Complex64>, Vec<Complex64>) {
    let n = f.len();

    let w = ROOTS_DICT.get(&(n as i32)).unwrap();

    let mut f0 = vec![Complex64::new(0.0, 0.0); n / 2];
    let mut f1 = vec![Complex64::new(0.0, 0.0); n / 2];

    for i in 0..(n / 2) {
        f0[i] = 0.5 * (f[2 * i] + f[2 * i + 1]);
        f1[i] = 0.5 * (f[2 * i] - f[2 * i + 1]) * w[2 * i].conj();
    }

    (f0, f1)
}

/*
    Compute the FFT of a polynomial mod (x^n + 1)

    Params:
        - f: A polynomial in coefficient representation

    Returns:
        - A polynomial in FFT representation
*/
fn fft(f: Vec<Complex64>) -> Vec<Complex64> {
    todo!();
}
