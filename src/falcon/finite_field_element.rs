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

use std::ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign};

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

impl Default for FiniteFieldElem {
    fn default() -> Self {
        FiniteFieldElem::new(0)
    }
}

impl Mul for FiniteFieldElem {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let new_val: u32 = (self.value as u64 * other.value as u64 % Q as u64) as u32;
        Self { value: new_val }
    }
}

impl Add for FiniteFieldElem {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let new_val: i64 = (self.value as i64 + other.value as i64) % Q as i64;
        Self {
            value: new_val as u32,
        }
    }
}

impl Sub for FiniteFieldElem {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let new_val: u32 = (self.value as i64 - other.value as i64 + Q as i64) as u32 % Q as u32;
        Self { value: new_val }
    }
}

impl AddAssign for FiniteFieldElem {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for FiniteFieldElem {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for FiniteFieldElem {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let neg_val = (Q as i32 - self.value as i32) % Q as i32;
        Self {
            value: neg_val as u32,
        }
    }
}

impl FiniteFieldElem {
    pub fn new(val: u32) -> Self {
        FiniteFieldElem {
            value: val % Q as u32,
        }
    }

    pub fn neg(&mut self) {
        let neg_val = (Q as i32 - self.value as i32) % Q as i32;

        self.value = neg_val as u32;
    }
}
