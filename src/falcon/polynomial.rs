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
use num_bigint::BigInt;

#[derive(Debug)]
struct Polynomial{
    coefficients: Vec<BigInt>
}

impl Polynomial{
   /*
   Splits a polynomial into two polynomials, one containing the even indices, and the other containing the odd.

   Used in the initial step of NTT computation to increase efficiency by allowing recursive calls to smaller
   degree polynomials.

   todo: may require padding of zeros on p to ensure that the number of coefficients is even.
    */
    fn split(p: &Polynomial) -> Vec<Polynomial>{
        todo!()
    }

    /*
    Merges two polynomials into a single polynomial by interleaving their coefficients.

    Used after NTT computations to recombine the polynomial. Called in reverse order to that of
    split(). merge() is also called recursively in this manner.
     */
    fn merge(v: Vec<Polynomial>) -> Polynomial{
        todo!()
    }
}