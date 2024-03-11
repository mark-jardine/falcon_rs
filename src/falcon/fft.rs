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

use super::polynomial::Polynomial;
use super::{fft_consts::ROOTS_DICT, finite_field_element::FiniteFieldElem};
use num_complex::{Complex, Complex64};
use std::vec;

#[derive(Debug)]
enum FftError {
    ZeroLengthPolynomial,
    InvalidRootIndex,
    FailedToSplitPolynomial,
    FailedToMergePolynomial, // InvNttFailed,
}

/*
    Split a polynomial f into two polynomials

    Params:
        - f: a polynomial in FFT representation

    Returns:
        - a Tuple, containing two polynomials in FFT representation

    Reference: Algorithm 1 splitfft(FFT(f)) in Falcon's specification.
*/
fn split_fft(f: Vec<Complex64>) -> Result<(Vec<Complex64>, Vec<Complex64>), FftError> {
    let length: usize = f.len();

    // Precomputed roots of unity associated with FFT of size `length`
    let w: &Vec<Complex<f64>> = match ROOTS_DICT.get(&length) {
        Some(val) => val,
        None => return Err(FftError::InvalidRootIndex),
    };

    let mut f0 = vec![Complex64::new(0.0, 0.0); length / 2];
    let mut f1 = vec![Complex64::new(0.0, 0.0); length / 2];

    for i in 0..(length / 2) {
        f0[i] = 0.5 * (f[2 * i] + f[2 * i + 1]);
        f1[i] = 0.5 * (f[2 * i] - f[2 * i + 1]) * w[2 * i].conj();
    }

    Ok((f0, f1))
}
/*
    Merge two or three polynomials into a single polynomial f.

    Params:
        f_vec_fft: a vector of polynomials in FFT representation

    Reference: algorithm 2 (mergefft_2) in Falcon's specification.
*/
fn merge_fft(f_vec_fft: Vec<Vec<Complex64>>) -> Result<Vec<Complex64>, FftError> {
    let f0_fft: Vec<Complex64> = f_vec_fft[0].clone();
    let f1_fft: Vec<Complex64> = f_vec_fft[1].clone();

    let length: usize = 2 * f0_fft.len();

    let w: &Vec<Complex<f64>> = match ROOTS_DICT.get(&length) {
        Some(val) => val,
        None => return Err(FftError::InvalidRootIndex),
    };

    let mut f_fft: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); length];

    for i in 0..(length / 2) {
        f_fft[2 * i] = f0_fft[i] + w[2 * i] * f1_fft[i];
        f_fft[2 * i + 1] = f0_fft[i] - w[2 * i] * f1_fft[i]
    }

    Ok(f_fft)
}

/*
    Compute the FFT of a polynomial mod (x^n + 1)

    Params:
        - f: A polynomial in coefficient representation

    Returns:
        - A polynomial in FFT representation
*/
fn fft(f: Polynomial) -> Result<Vec<Complex64>, FftError> {
    let length: usize = f.coefficients.len();
    let mut f_fft: Vec<Complex64> = vec![];

    //Get u32 values as Vec<u32> from Polyomial's Vec<FiniteFieldElem> `coefficients`
    let f_coeff: Vec<u32> = f.coefficients.iter().map(|x| x.value).collect();

    if length > 2 {
        let (f0, f1) = Polynomial::split(&f);
        let f0_fft = fft(f0)?;
        let f1_fft = fft(f1)?;
        f_fft = match merge_fft(vec![f0_fft, f1_fft]) {
            Ok(val) => val,
            Err(_e) => return Err(FftError::FailedToMergePolynomial),
        }
    } else if length == 2 {
        f_fft = vec![Complex64::new(0.0, 0.0); length];
        f_fft[0] = Complex64::new(f_coeff[0] as f64, 1.0 * f_coeff[1] as f64);
        f_fft[1] = Complex64::new(f_coeff[0] as f64, -1.0 * f_coeff[1] as f64);
    }

    Ok(f_fft)
}

fn inv_fft(f_fft: Vec<Complex64>) -> Result<Vec<f64>, FftError> {
    let length: usize = f_fft.len();

    let mut f: Vec<f64> = vec![];

    if length > 2 {
        let (f0_fft, f1_fft) = split_fft(f_fft)?;
        let f0 = inv_fft(f0_fft)?;
        let f1 = inv_fft(f1_fft)?;
        f = Polynomial::merge_fp(f0, f1);
    } else if length == 2 {
        f = vec![0.0; length];
        f[0] = f_fft[0].re;
    }

    Ok(f)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::falcon::finite_field_element::FiniteFieldElem;

    #[test]
    fn test_fft_degree_preservation() {
        let f = Polynomial::new(vec![
            FiniteFieldElem::new(1),
            FiniteFieldElem::new(2),
            FiniteFieldElem::new(3),
            FiniteFieldElem::new(4),
        ]);

        let f_fft = fft(f).expect("fft() failed in test_fft_degree_preservation().");

        // Assuming fft() pads the input to the nearest power of 2
        let expected_length = 4;
        assert_eq!(
            f_fft.len(),
            expected_length,
            "FFT should preserve the degree after padding."
        );
    }
}
