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
/*
    A polynomial that accepts a vector of generic T.

    For coefficient format, `FiniteFieldElem` is used
    For FFT format, `Complex64` is used.
*/

use super::{
    fft::{fft, inv_fft},
    finite_field_element::FiniteFieldElem,
    ntrugen::{self, bitsize},
};
use num_bigint::BigInt;
use num_complex::Complex64;
use std::{
    error::Error,
    ops::{Add, Div, Mul, Sub},
};

#[derive(Debug, Clone)]
pub struct Polynomial<T> {
    pub coefficients: Vec<T>,
}

impl<T> Polynomial<T>
where
    T: Clone,
{
    pub fn new(_coefficients: Vec<T>) -> Self {
        Polynomial {
            coefficients: _coefficients,
        }
    }
}

impl<T> FromIterator<T> for Polynomial<T>
where
    T: Clone + Copy + Default + Mul + Add + Sub,
{
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let coefficients = iter.into_iter().collect::<Vec<_>>();
        Polynomial::<T>::new(coefficients)
    }
}

impl Polynomial<FiniteFieldElem> {
    /*
       Galois conjugate of an element a in Q[x] / (x ** n + 1).
       Here, the Galois conjugate of a(x) is simply a(-x).

       Params:
            - a reference to self: a Polynomial<FiniteFieldElem>

        Returns:
            - a Polynomial<FiniteFieldElem> representing the negated polynomial
    */
    pub fn galois_conjugate(&self) -> Self {
        let conjugate_coeffs: Vec<FiniteFieldElem> = self
            .coefficients
            .iter()
            .enumerate()
            .map(|(i, &coeff)| if i % 2 == 0 { coeff } else { -coeff })
            .collect();

        Polynomial::new(conjugate_coeffs)
    }

    /*
       Project an element a of Q[x] / (x ** n + 1) onto Q[x] / (x ** (n // 2) + 1).
       n must be a power-of-two.

      Params:
           - a reference to self: a Polynomial<FiniteFieldElem>

       Returns:
           - a Polynomial<FiniteFieldElem>
    */
    pub fn field_norm(&self) -> Self {
        let n = self.coefficients.len() / 2;
        let mut even_coeffs: Vec<FiniteFieldElem> = vec![FiniteFieldElem::default(); n];
        let mut odd_coeffs: Vec<FiniteFieldElem> = vec![FiniteFieldElem::default(); n];

        for (i, &coeff) in self.coefficients.iter().enumerate() {
            if i % 2 == 0 {
                even_coeffs.push(coeff);
            } else {
                odd_coeffs.push(coeff);
            }
        }
        let even_coeffs_sq = karamul(&even_coeffs, &even_coeffs);
        let odd_coeffs_sq = karamul(&odd_coeffs, &odd_coeffs);

        let mut res = even_coeffs_sq.clone();
        res.iter_mut()
            .skip(1)
            .enumerate()
            .for_each(|(index, curr)| *curr -= even_coeffs_sq[index - 1]);
        res[0] = odd_coeffs_sq[n - 1];

        Polynomial::new(res)
    }

    /*
       Lift an element a of Q[x] / (x ** (n // 2) + 1) up to Q[x] / (x ** n + 1).
       The lift of a(x) is simply a(x ** 2) seen as an element of Q[x] / (x ** n + 1).

       Params:
          - a reference to self: a Polynomial<FiniteFieldElem>

       Returns:
          - a Polynomial<FiniteFieldElem> with double the degree
    */
    fn lift(&self) -> Self {
        let mut res = vec![FiniteFieldElem::default(); self.coefficients.len() * 2];

        self.coefficients
            .iter()
            .enumerate()
            .for_each(|(index, coeff)| res[2 * index] = *coeff);

        Polynomial::new(res)
    }

    /*
        Splits a polynomial into two polynomials, one containing the even indices, and the other containing the odd.

        Used in the initial step of NTT computation to increase efficiency by allowing recursive calls to smaller
        degree polynomials.

        This will produce incorrect values in the case that the number of coefficients
        in p is not even.
    */
    pub fn split(
        p: &Polynomial<FiniteFieldElem>,
    ) -> (Polynomial<FiniteFieldElem>, Polynomial<FiniteFieldElem>) {
        let length: usize = p.coefficients.len();

        let mut f0_coeffs: Vec<FiniteFieldElem> = Vec::new();
        let mut f1_coeffs: Vec<FiniteFieldElem> = Vec::new();

        for i in 0..(length / 2) {
            f0_coeffs.push(p.coefficients[2 * i]);
            f1_coeffs.push(p.coefficients[2 * i + 1]);
        }

        let f0: Polynomial<FiniteFieldElem> = Polynomial::new(f0_coeffs);
        let f1: Polynomial<FiniteFieldElem> = Polynomial::new(f1_coeffs);

        (f0, f1)
    }
    /*
       Merges two polynomials into a single polynomial by interleaving their coefficients.

       Used after NTT computations to recombine the polynomial. Called in reverse order to that of
       split(). merge() is also called recursively in this manner.
    */
    pub fn merge(f_vec: Vec<Polynomial<FiniteFieldElem>>) -> Polynomial<FiniteFieldElem> {
        let f0 = f_vec.first().unwrap();
        let f1 = f_vec.get(1).unwrap();
        let length: usize = f0.coefficients.len() * 2;
        let mut f: Polynomial<FiniteFieldElem> =
            Polynomial::new(vec![FiniteFieldElem::default(); length]);

        for i in 0..length / 2 {
            f.coefficients[2 * i] = f0.coefficients[i];
            f.coefficients[2 * i + 1] = f1.coefficients[i];
        }

        f
    }
}

impl Polynomial<Complex64> {
    /*
        Splits a fft polynomial into two polynomials, one containing the even indices, and the other containing the odd.

        Floating point version of Polynomial::split(), used for FFT.
    */
    pub fn split_fp(f: &Polynomial<Complex64>) -> (Polynomial<Complex64>, Polynomial<Complex64>) {
        let mut f0: Polynomial<Complex64> = Polynomial::new(vec![]);
        let mut f1: Polynomial<Complex64> = Polynomial::new(vec![]);

        for (i, &value) in f.coefficients.iter().enumerate() {
            if i % 2 == 0 {
                f0.coefficients.push(value);
            } else {
                f1.coefficients.push(value);
            }
        }

        (f0, f1)
    }

    pub fn adjoint(&mut self) -> Self {
        Polynomial::new(self.coefficients.iter().map(|c| c.conj()).collect())
    }
}

impl Add for Polynomial<Complex64> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Polynomial::new(
            self.coefficients
                .iter()
                .zip(rhs.coefficients)
                .map(|(s, rhs)| s + rhs)
                .collect(),
        )
    }
}
impl Sub for Polynomial<Complex64> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Polynomial::new(
            self.coefficients
                .iter()
                .zip(rhs.coefficients)
                .map(|(s, rhs)| s - rhs)
                .collect(),
        )
    }
}

impl Mul for Polynomial<Complex64> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Polynomial::new(
            self.coefficients
                .iter()
                .zip(rhs.coefficients)
                .map(|(s, rhs)| s * rhs)
                .collect(),
        )
    }
}

impl Div for Polynomial<Complex64> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Polynomial::new(
            self.coefficients
                .iter()
                .zip(rhs.coefficients)
                .map(|(s, rhs)| s / rhs)
                .collect(),
        )
    }
}

impl Polynomial<f64> {
    /*
       Merges two fft polynomials into a single polynomial by interleaving their coefficients.

       Floating point version of Polynomial::merge(), used for FFT.
    */
    pub fn merge_fp(f0: Polynomial<f64>, f1: Polynomial<f64>) -> Polynomial<f64> {
        let length = f0.coefficients.len() * 2;
        let mut f_merged: Polynomial<f64> = Polynomial::new(vec![f64::default(); length]);

        for i in 0..f0.coefficients.len() {
            f_merged.coefficients[2 * i] = f0.coefficients[i];
            f_merged.coefficients[2 * i + 1] = f1.coefficients[i];
        }

        f_merged
    }
}

/*
    Karatsuba multiplication algorithm for Polynomials

    Params:
        - a, slice of coefficients
        - b, a slice of coefficients

    Returns:
        - a vector of coefficients
*/
pub fn karatsuba<T>(a: &[T], b: &[T]) -> Vec<T>
where
    T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Default + Clone,
{
    let n = a.len();
    if n <= 1 {
        return vec![a[0].clone() * b[0].clone()];
    }

    let mid = n / 2;
    let (a0, a1) = a.split_at(mid);
    let (b0, b1) = b.split_at(mid);

    let a0b0 = karatsuba(a0, b0);
    let a1b1 = karatsuba(a1, b1);

    let a0_plus_a1: Vec<T> = a0
        .iter()
        .cloned()
        .zip(a1.iter().cloned())
        .map(|(x, y)| x + y)
        .collect();
    let b0_plus_b1: Vec<T> = b0
        .iter()
        .cloned()
        .zip(b1.iter().cloned())
        .map(|(x, y)| x + y)
        .collect();

    let middle = karatsuba(&a0_plus_a1, &b0_plus_b1)
        .into_iter()
        .zip(a0b0.iter().cloned())
        .zip(a1b1.iter().cloned())
        .map(|((x, y), z)| x - y - z)
        .collect::<Vec<T>>();

    // Combine the results: a0b0 + middle + a1b1
    let mut result = vec![T::default(); n * 2];
    for i in 0..n {
        result[i] = result[i].clone() + a0b0.get(i).unwrap_or(&T::default()).clone();
        result[i + mid] = result[i + mid].clone() + middle.get(i).unwrap_or(&T::default()).clone();
        result[i + n] = result[i + n].clone() + a1b1.get(i).unwrap_or(&T::default()).clone();
    }

    result
}

/*
    Karatsuba multiplication, followed by reduction mod (x ** n + 1).

    Params:
        - a, slice of coefficients
        - b, a slice of coefficients

    Returns:
        - a vector of coefficients
        -
*/
pub fn karamul<T>(a: &[T], b: &[T]) -> Vec<T>
where
    T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Default + Clone,
{
    let ab = karatsuba(a, b);
    let mut abr: Vec<T> = vec![T::default(); ab.len() / 2];
    let length = a.len();

    for i in 0..length {
        abr[i] = ab[i].clone() - ab[i + length].clone();
    }

    abr
}

impl Polynomial<BigInt> {
    /*
       Reduce (F, G) relatively to (f, g).

       This is done via Babai's reduction.
       (F, G) <-- (F, G) - k * (f, g), where k = round((F f* + G g*) / (f f* + g g*)).
       Corresponds to algorithm 7 (Reduce) of Falcon's documentation.

       This function modifies F, G in place

       Params:
            - f, g Polynomial<BigInt>
            - F, G &mut Polynomial<BigInt>

        Returns:
            - Either void(meaning F, G have been modified in place), or an Error

    */
    pub fn reduce(
        f: Polynomial<BigInt>,
        g: Polynomial<BigInt>,
        capf: &mut Polynomial<BigInt>,
        capg: &mut Polynomial<BigInt>,
    ) -> Result<(), Box<dyn Error>> {
        let len = f.coefficients.len();
        let f_min = f.coefficients.iter().min().unwrap();
        let f_max = f.coefficients.iter().max().unwrap();
        let g_min = g.coefficients.iter().min().unwrap();
        let g_max = g.coefficients.iter().max().unwrap();

        let size = vec![
            53u64,
            ntrugen::bitsize(f_min.to_owned()),
            ntrugen::bitsize(f_max.to_owned()),
            ntrugen::bitsize(g_min.to_owned()),
            ntrugen::bitsize(g_max.to_owned()),
        ]
        .into_iter()
        .max()
        .ok_or("Could not find size in reduce() in polynomial.rs")?;

        let f_adjust_coeffs: Vec<Complex64> = f
            .coefficients
            .iter()
            .map(|elt| elt >> (size - 53))
            .map(|elt| Complex64::new(i64::try_from(elt).unwrap() as f64, 0.0))
            .collect();

        let g_adjust_coeffs: Vec<Complex64> = g
            .coefficients
            .iter()
            .map(|elt| elt >> (size - 53))
            .map(|elt| Complex64::new(i64::try_from(elt).unwrap() as f64, 0.0))
            .collect();

        let f_adjust = Polynomial::new(f_adjust_coeffs);
        let g_adjust = Polynomial::new(g_adjust_coeffs);
        let mut f_adj_fft = fft(f_adjust)?;
        let mut g_adj_fft = fft(g_adjust)?;

        loop {
            let capf_min = capf.coefficients.iter().min().unwrap();
            let capf_max = capf.coefficients.iter().max().unwrap();
            let capg_min = capg.coefficients.iter().min().unwrap();
            let capg_max = capf.coefficients.iter().max().unwrap();

            let cap_size = vec![
                53u64,
                bitsize(capf_min.to_owned()),
                bitsize(capf_max.to_owned()),
                bitsize(capg_min.to_owned()),
                bitsize(capg_max.to_owned()),
            ]
            .into_iter()
            .max()
            .ok_or("Could not find Size in reduce() in polynomial.rs")?;

            if cap_size < size {
                break;
            }

            let capf_adjust_coeffs: Vec<Complex64> = capf
                .coefficients
                .iter()
                .map(|elt| elt >> (size - 53))
                .map(|elt| Complex64::new(i64::try_from(elt).unwrap() as f64, 0.0))
                .collect();

            let capg_adjust_coeffs: Vec<Complex64> = capg
                .coefficients
                .iter()
                .map(|elt| elt >> (size - 53))
                .map(|elt| Complex64::new(i64::try_from(elt).unwrap() as f64, 0.0))
                .collect();

            let capf_adjust = Polynomial::new(capf_adjust_coeffs);
            let capg_adjust = Polynomial::new(capg_adjust_coeffs);
            let capf_adj_fft = fft(capf_adjust)?;
            let capg_adj_fft = fft(capg_adjust)?;

            let den_fft = (f_adj_fft.clone() * f_adj_fft.adjoint())
                + (g_adj_fft.clone() * g_adj_fft.adjoint());
            let num_fft =
                (capf_adj_fft * f_adj_fft.adjoint()) + (capg_adj_fft * g_adj_fft.adjoint());
            let k_fft = num_fft / den_fft;
            let k_fp = inv_fft(k_fft)?;
            let k: Polynomial<i64> = k_fp.coefficients.iter().map(|f| f.round() as i64).collect();

            let k_all_zeroes: bool = k
                .coefficients
                .iter()
                .filter(|e| **e == 0)
                .collect::<Vec<&i64>>()
                .is_empty();
            if k_all_zeroes {
                break;
            }
            let k_bigint_coeffs: Vec<BigInt> = k
                .coefficients
                .iter()
                .map(|&num| BigInt::from(num))
                .collect();

            // TODO: Following can be replaced with Toom-Cook when optmising
            let fk = karamul(&f.coefficients, &k_bigint_coeffs);
            let gk = karamul(&g.coefficients, &k_bigint_coeffs);

            for i in 0..len {
                capf.coefficients[i] -= &fk[i] << (cap_size - size);
                capg.coefficients[i] -= &gk[i] << (cap_size - size);
            }
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn test_karatsuba() {
        // Test case for FiniteFieldElem
        {
            let a_coeffs = vec![FiniteFieldElem::new(1), FiniteFieldElem::new(2)];
            let b_coeffs = vec![FiniteFieldElem::new(3), FiniteFieldElem::new(4)];
            let a = Polynomial::new(a_coeffs);
            let b = Polynomial::new(b_coeffs);

            let result = karatsuba(&a.coefficients, &b.coefficients);
            let expected_coeffs = vec![
                FiniteFieldElem::new(3),
                FiniteFieldElem::new(10),
                FiniteFieldElem::new(8),
                FiniteFieldElem::new(0), // karatsuba() will pad with 0 on odd-length polynomials
            ];
            let expected = Polynomial::new(expected_coeffs);

            assert_eq!(result, expected.coefficients);
        }

        // Test case for Complex64
        {
            let a_coeffs = vec![Complex64::new(1.0, -1.0), Complex64::new(2.0, 3.0)];
            let b_coeffs = vec![Complex64::new(4.0, 2.0), Complex64::new(-1.0, 1.0)];
            let result = karatsuba(&a_coeffs, &b_coeffs);

            let expected_coeffs = vec![
                Complex64::new(6.0, -2.0),
                Complex64::new(2.0, 18.0),
                Complex64::new(-5.0, -1.0),
            ];
            let expected_poly = Polynomial::new(expected_coeffs);

            assert!(result
                .iter()
                .zip(expected_poly.coefficients.iter())
                .all(|(a, b)| (a - b).norm() < 1e-9));
        }
    }
}
