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
fn fft(f: Vec<Complex64>) -> Result<Vec<Complex64>, FftError> {
    let length: usize = f.len();
    let mut f_fft: Vec<Complex64> = vec![];
    println!("f:{f:?}\n");

    if length > 2 {
        let (f0, f1) = Polynomial::split_fp(&f);
        let f0_fft = fft(f0)?;
        let f1_fft = fft(f1)?;
        f_fft = merge_fft(vec![f0_fft, f1_fft])?;
    } else if length == 2 {
        // Base case handling as before
        f_fft = vec![Complex64::new(0.0, 0.0); length];
        f_fft[0] = f[0] + Complex64::new(0.0, 1.0) * f[1];
        f_fft[1] = f[0] - Complex64::new(0.0, 1.0) * f[1];
    }

    println!("\nf_fft:{f_fft:?}");
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

    #[test]
    fn test_fft_degree_preservation() {
        let f = vec![
            Complex64::new(1.0, 1.0),
            Complex64::new(2.0, 2.0),
            Complex64::new(3.0, 3.0),
            Complex64::new(4.0, 4.0),
        ];

        let f_fft = fft(f).expect("fft() failed in test_fft_degree_preservation().");

        // Assuming fft() pads the input to the nearest power of 2
        let expected_length = 4;
        assert_eq!(
            f_fft.len(),
            expected_length,
            "FFT should preserve the degree after padding."
        );
    }

    #[test]
    fn test_fft_basic() {
        let inputs_and_expected = vec![
            (
                vec![
                    Complex64::new(1.0, 0.0),
                    Complex64::new(2.0, 0.0),
                    Complex64::new(3.0, 0.0),
                    Complex64::new(4.0, 0.0),
                ],
                vec![
                    Complex64::new(-0.41421356237309204, 7.242640687119286),
                    Complex64::new(2.4142135623730923, -1.2426406871192857),
                    Complex64::new(-0.41421356237309204, -7.242640687119286),
                    Complex64::new(2.4142135623730923, 1.2426406871192857),
                ],
            ),
            (
                vec![
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                ],
                vec![
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                ],
            ),
        ];

        for (input, expected) in inputs_and_expected {
            let result = fft(input); // Assuming your Rust FFT function is named `fft` and returns a `Result<Vec<Complex64>, Error>`
            assert_eq!(result.unwrap(), expected);
        }
    }

    fn assert_complex_vec_eq(result: &[Complex64], expected: &[Complex64], epsilon: f64) {
        assert_eq!(result.len(), expected.len(), "Vector lengths differ");
        for (i, (r, e)) in result.iter().zip(expected.iter()).enumerate() {
            assert!(
                (r.re - e.re).abs() < epsilon && (r.im - e.im).abs() < epsilon,
                "Mismatch at index {}: expected {:?}, got {:?}",
                i,
                e,
                r
            );
        }
    }

    #[test]
    fn test_fft_zero_polynomial() {
        let input = vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ];
        let expected_output = vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ];
        let result = fft(input).unwrap();
        assert_complex_vec_eq(&result, &expected_output, 1e-9);
    }

    #[test]
    fn test_fft_single_non_zero_coefficient() {
        let input = vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(10.0, 0.0),
        ];
        let expected_output = vec![
            Complex64::new(-7.071067811865475, 7.071067811865475),
            Complex64::new(7.071067811865475, -7.071067811865475),
            Complex64::new(-7.071067811865475, -7.071067811865475),
            Complex64::new(7.071067811865475, 7.071067811865475),
        ];
        let result = fft(input).unwrap();
        assert_complex_vec_eq(&result, &expected_output, 1e-9);
    }

    #[test]
    fn test_fft_all_ones_polynomial() {
        let input = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        let expected_output = vec![
            Complex64::new(1.0, 2.414213562373095),
            Complex64::new(1.0, -0.4142135623730949),
            Complex64::new(1.0, -2.414213562373095),
            Complex64::new(1.0, 0.4142135623730949),
        ];
        let result = fft(input).unwrap();
        assert_complex_vec_eq(&result, &expected_output, 1e-9);
    }

    #[test]
    fn test_fft_alternating_sign_polynomial() {
        let input = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(-1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(-1.0, 0.0),
        ];
        let expected_output = vec![
            Complex64::new(1.0, -0.4142135623730949),
            Complex64::new(1.0, 2.414213562373095),
            Complex64::new(1.0, 0.4142135623730949),
            Complex64::new(1.0, -2.414213562373095),
        ];

        let result = fft(input).unwrap();
        assert_complex_vec_eq(&result, &expected_output, 1e-9);
    }
}
