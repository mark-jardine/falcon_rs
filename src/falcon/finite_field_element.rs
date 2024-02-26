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
    Q is the modulus used in Falcon (12 * 1024 + 1)
*/
pub const Q: u16 = 12289;

/*
   32-bit unsigned integer wrapper used to represent a finite field element.
   Handles the finite field arithmetic.
*/
#[derive(Debug, Clone, Copy)]
pub struct FiniteFieldElem {
    pub value: u32,
}

impl PartialEq for FiniteFieldElem {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

//todo: maybe change val to be u32, likewise for Q
impl FiniteFieldElem {
    pub fn new(val: u32) -> Self {
        FiniteFieldElem {
            value: val % Q as u32,
        }
    }

    pub fn mult(&mut self, other: &FiniteFieldElem) {
        let new_val: i64 = (self.value as i64 * other.value as i64) % Q as i64;

        self.value = new_val as u32;
    }

    pub fn add(&mut self, other: &FiniteFieldElem) {
        let new_val: i64 = (self.value as i64 + other.value as i64) % Q as i64;

        self.value = new_val as u32;
    }

    pub fn sub(&mut self, other: &FiniteFieldElem) {
        self.value = (self.value as i64 - other.value as i64 + Q as i64) as u32 % Q as u32;
    }

    pub fn neg(&mut self) {
        let neg_val = (Q as i32 - self.value as i32) % Q as i32;

        self.value = neg_val as u32;
    }
}
