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
use num_complex::{Complex, Complex64};
use std::vec;

#[derive(Debug)]
pub enum FftError {
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
fn split_fft(
    f: Polynomial<Complex64>,
) -> Result<(Polynomial<Complex64>, Polynomial<Complex64>), FftError> {
    let length: usize = f.coefficients.len();

    // Precomputed roots of unity associated with FFT of size `length`
    let w: &Vec<Complex<f64>> = match ROOTS_DICT.get(&length) {
        Some(val) => val,
        None => return Err(FftError::InvalidRootIndex),
    };

    let mut f0 = Polynomial::new(vec![Complex64::new(0.0, 0.0); length / 2]);
    let mut f1 = Polynomial::new(vec![Complex64::new(0.0, 0.0); length / 2]);

    for i in 0..(length / 2) {
        f0.coefficients[i] = 0.5 * (f.coefficients[2 * i] + f.coefficients[2 * i + 1]);
        f1.coefficients[i] =
            0.5 * (f.coefficients[2 * i] - f.coefficients[2 * i + 1]) * w[2 * i].conj();
    }

    Ok((f0, f1))
}
/*
    Merge two or three polynomials into a single polynomial f.

    Params:
        f_vec_fft: a vector of polynomials in FFT representation

    Reference: algorithm 2 (mergefft_2) in Falcon's specification.
*/
fn merge_fft(f_vec_fft: Vec<Polynomial<Complex64>>) -> Result<Polynomial<Complex64>, FftError> {
    let f0_fft: &Polynomial<Complex<f64>> = &f_vec_fft[0];
    let f1_fft: &Polynomial<Complex<f64>> = &f_vec_fft[1];

    let length: usize = 2 * f0_fft.coefficients.len();

    let w: &Vec<Complex<f64>> = match ROOTS_DICT.get(&length) {
        Some(val) => val,
        None => return Err(FftError::InvalidRootIndex),
    };

    let mut f_fft: Polynomial<Complex<f64>> =
        Polynomial::new(vec![Complex64::new(0.0, 0.0); length]);

    for i in 0..(length / 2) {
        f_fft.coefficients[2 * i] = f0_fft.coefficients[i] + w[2 * i] * f1_fft.coefficients[i];
        f_fft.coefficients[2 * i + 1] = f0_fft.coefficients[i] - w[2 * i] * f1_fft.coefficients[i];
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
pub fn fft(f: Polynomial<Complex64>) -> Result<Polynomial<Complex64>, FftError> {
    let length: usize = f.coefficients.len();
    let mut f_fft: Polynomial<Complex64> = Polynomial::new(vec![]);

    if length > 2 {
        let (f0, f1) = Polynomial::<Complex64>::split_fp(&f);
        let f0_fft = fft(f0)?;
        let f1_fft = fft(f1)?;
        f_fft = merge_fft(vec![f0_fft, f1_fft])?;
    } else if length == 2 {
        // Base case handling as before
        f_fft = Polynomial::new(vec![Complex64::new(0.0, 0.0); length]);
        f_fft.coefficients[0] = f.coefficients[0] + Complex64::new(0.0, 1.0) * f.coefficients[1];
        f_fft.coefficients[1] = f.coefficients[0] - Complex64::new(0.0, 1.0) * f.coefficients[1];
    }

    println!("\nf_fft:{f_fft:?}");
    Ok(f_fft)
}

fn inv_fft(f_fft: Polynomial<Complex64>) -> Result<Polynomial<f64>, FftError> {
    let length: usize = f_fft.coefficients.len();

    let mut f: Polynomial<f64> = Polynomial::new(vec![]);

    if length > 2 {
        let (f0_fft, f1_fft) = split_fft(f_fft)?;
        let f0 = inv_fft(f0_fft)?;
        let f1 = inv_fft(f1_fft)?;
        f = Polynomial::<f64>::merge_fp(f0, f1);
    } else if length == 2 {
        f.coefficients = vec![0.0; length];
        f.coefficients[0] = f_fft.coefficients[0].re;
    }

    Ok(f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fft_degree_preservation() {
        let f = Polynomial::new(vec![
            Complex64::new(1.0, 1.0),
            Complex64::new(2.0, 2.0),
            Complex64::new(3.0, 3.0),
            Complex64::new(4.0, 4.0),
        ]);

        let f_fft = fft(f).expect("fft() failed in test_fft_degree_preservation().");

        // Assuming fft() pads the input to the nearest power of 2
        let expected_length = 4;
        assert_eq!(
            f_fft.coefficients.len(),
            expected_length,
            "FFT should preserve the degree after padding."
        );
    }

    #[test]
    fn test_fft_basic() {
        let inputs_and_expected = vec![
            (
                Polynomial::new(vec![
                    Complex64::new(1.0, 0.0),
                    Complex64::new(2.0, 0.0),
                    Complex64::new(3.0, 0.0),
                    Complex64::new(4.0, 0.0),
                ]),
                Polynomial::new(vec![
                    Complex64::new(-0.41421356237309204, 7.242640687119286),
                    Complex64::new(2.4142135623730923, -1.2426406871192857),
                    Complex64::new(-0.41421356237309204, -7.242640687119286),
                    Complex64::new(2.4142135623730923, 1.2426406871192857),
                ]),
            ),
            (
                Polynomial::new(vec![
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                ]),
                Polynomial::new(vec![
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                ]),
            ),
        ];

        for (input, expected) in inputs_and_expected {
            let result = fft(input);
            assert_eq!(result.unwrap().coefficients, expected.coefficients);
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
        let input = Polynomial::new(vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]);
        let expected_output = Polynomial::new(vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]);
        let result = fft(input).unwrap();
        assert_complex_vec_eq(&result.coefficients, &expected_output.coefficients, 1e-9);
    }

    #[test]
    fn test_fft_single_non_zero_coefficient() {
        let input = Polynomial::new(vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(10.0, 0.0),
        ]);
        let expected_output = Polynomial::new(vec![
            Complex64::new(-7.071067811865475, 7.071067811865475),
            Complex64::new(7.071067811865475, -7.071067811865475),
            Complex64::new(-7.071067811865475, -7.071067811865475),
            Complex64::new(7.071067811865475, 7.071067811865475),
        ]);
        let result = fft(input).unwrap();
        assert_complex_vec_eq(&result.coefficients, &expected_output.coefficients, 1e-9);
    }

    #[test]
    fn test_fft_all_ones_polynomial() {
        let input = Polynomial::new(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
        ]);
        let expected_output = Polynomial::new(vec![
            Complex64::new(1.0, 2.414213562373095),
            Complex64::new(1.0, -0.4142135623730949),
            Complex64::new(1.0, -2.414213562373095),
            Complex64::new(1.0, 0.4142135623730949),
        ]);
        let result = fft(input).unwrap();
        assert_complex_vec_eq(&result.coefficients, &expected_output.coefficients, 1e-9);
    }

    #[test]
    fn test_fft_alternating_sign_polynomial() {
        let input = Polynomial::new(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(-1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(-1.0, 0.0),
        ]);
        let expected_output = Polynomial::new(vec![
            Complex64::new(1.0, -0.4142135623730949),
            Complex64::new(1.0, 2.414213562373095),
            Complex64::new(1.0, 0.4142135623730949),
            Complex64::new(1.0, -2.414213562373095),
        ]);

        let result = fft(input).unwrap();
        assert_complex_vec_eq(&result.coefficients, &expected_output.coefficients, 1e-9);
    }
}
