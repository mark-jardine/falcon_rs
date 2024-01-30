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
struct FiniteFieldElem {
    value: u32
}

impl FiniteFieldElem{
    pub fn new(val: i32) -> Self{
        // Get sign, which will be either 1 or 0.
        // 1 means that input value >= 0
        // 0 means that input value < 0
        let gte_zero: bool = val >= 0;
        let gte_zero_int: i32 = gte_zero as i32;
        let sign: i32 = gte_zero_int - ((!gte_zero) as i32);

        // Apply modulus, if the value is negative, val_mod = 0
        let val_mod: i32 = sign * (sign * (val % Q as i32));

        // If the input value is negative, return Q, otherwise return the modulus-reduced value
        let f_elem: u32 = (val_mod as u32 + (Q as u32 * (1 - gte_zero_int as u32)));

        FiniteFieldElem{value: f_elem}
    }

    pub fn mult(&self, other: FiniteFieldElem) -> FiniteFieldElem{
        let new_val: u32 = (self.value * other.value) % Q as u32;

        FiniteFieldElem{value: new_val}
    }

    pub fn add(&self, other: FiniteFieldElem) -> FiniteFieldElem{
        let new_val: u32 = (self.value + other.value) % Q as u32;

        FiniteFieldElem{value: new_val}
    }

    pub fn sub(&self, other: FiniteFieldElem) -> FiniteFieldElem{
        let new_val: u32 = (self.value - other.value);

        // If the new value is negative, return Q, otherwise return the modulus-reduced value
        let gte_zero: u32 = (new_val >= 0) as u32;
        let mod_val: u32 = new_val + (Q as u32 * (1 - gte_zero) );

        FiniteFieldElem{value: mod_val}
    }

}