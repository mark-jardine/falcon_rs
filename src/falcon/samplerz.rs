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

use std::{os, ptr::null};

// Use rand::OsRng for high-quality randomness
use rand::{rngs::OsRng, Rng, RngCore};

// Upper bound on all the values of sigma
const MAX_SIGMA: f64 = 1.8205;
// INV_2SIGMA2 = 1 / (2 * (MAX_SIGMA ** 2)), where (MAX_SIGMA ** 2) has been pre-computed to be 3.31422025
const INV_2SIGMA2: f64 = 1.0 / (2.0 * 3.31422025);

// Precision of RCDT
const RCDT_PREC: u8 = 72;

// ln(2) and 1 / ln(2), with ln the natural logarithm
const LN2: f64 = 0.69314718056;
const ILN2: f64 = 1.44269504089;

/*
    RCDT is the reverse cumulative distribution table of a distribution that
    is very close to a half-Gaussian of parameter MAX_SIGMA.
*/
const RCDT: [u128; 18] = [
    3024686241123004913666,
    1564742784480091954050,
    636254429462080897535,
    199560484645026482916,
    47667343854657281903,
    8595902006365044063,
    1163297957344668388,
    117656387352093658,
    8867391802663976,
    496969357462633,
    20680885154299,
    638331848991,
    14602316184,
    247426747,
    3104126,
    28824,
    198,
    1,
];

/*
    C contains the coefficients of a polynomial that approximates exp(-x)
    More precisely, the value:
    (2 ** -63) * sum(C[12 - i] * (x ** i) for i in range(i))
    Should be very close to exp(-x).
    This polynomial is lifted from FACCT: https://doi.org/10.1109/TC.2019.2940949
*/
const C: [u64; 13] = [
    0x00000004741183A3,
    0x00000036548CFC06,
    0x0000024FDCBF140A,
    0x0000171D939DE045,
    0x0000D00CF58F6F84,
    0x000680681CF796E3,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    0x8000000000000000,
];

/*
    Sample z0 in {0, 1, ..., 18} with a distribution very close to the half-Gaussian D_{Z+, 0, MAX_SIGMA}.
*/
fn base_sampler<R: RngCore>(rng: &mut R) -> u8 {
    let random_bytes_needed = (RCDT_PREC + 7) / 8;
    let mut buf = [0u8; 16];

    rng.fill_bytes(&mut buf);
    let mut random_value: u128 = 0;
    for i in 0..random_bytes_needed {
        random_value |= (buf[i as usize] as u128) << (8 * i);
    }

    let mut z0: u8 = 0;
    for &elt in RCDT.iter() {
        if random_value < elt {
            z0 += 1;
        } else {
            break;
        }
    }

    z0
}

/*
    Compute an approximation of 2^63 * ccs * exp(-x).

    Params:
    - a floating-point number x
    - a scaling factor ccs
    x > 0, ccs > 0

    Returns:
    - an integral approximation of 2^63 * ccs * exp(-x).
*/

fn approxexp(x: f64, ccs: f64) -> u64 {
    assert!(
        x > 0.0 && ccs > 0.0,
        "approxexp(): Both x and ccs must be positive.\nx:{x}, ccs:{ccs}"
    );

    // Convert x to a 64-bit integer representation for scaling
    let two_pow_63 = 1u64 << 63;
    let z: u64 = (x * two_pow_63 as f64) as u64;
    // Start with the first coefficient for polynomial evaluation
    let mut y = C[0];
    // Polynomial approximation calculation
    for &elt in &C[1..] {
        // Horner's Method
        let z_mul_y = z as u128 * y as u128;
        y = elt - (z_mul_y >> 63) as u64;
    }

    // Apply the scaling factor ccs
    let z_scaled: u128 = (f64::floor(ccs * (two_pow_63 as f64))) as u128;

    (z_scaled * y as u128 >> 63) as u64
}

/*
    Return a single bit, equal to 1 with probability ~ ccs * exp(-x).
    x > 0, ccs > 0
*/
fn berexp<R: RngCore>(x: f64, ccs: f64, rng: &mut R) -> bool {
    assert!(
        x > 0.0 && ccs > 0.0,
        "berexp(): Both x and ccs must be positive.\nx:{x}, ccs:{ccs}"
    );

    let mut s: i64 = (x * ILN2) as i64;
    let r: f64 = x - s as f64 * LN2;
    s = i64::min(s, 63);

    // z = (approxexp(r, ccs) - 1) >> s
    let z = (((approxexp(r, ccs) as u128) - 1) >> s) as u64;

    let mut w: i32 = 0;
    for i in (0..=56).rev().step_by(8) {
        let p: u8 = rng.gen::<u8>();

        // w = p - ((z >> i) & 0xFF)
        let last_8_bits = (z >> i) & 0xFF;
        w = p as i32 - last_8_bits as i32;

        if w != 0 {
            break;
        }
    }

    w < 0
}

/*
    Given floating-point values mu, sigma (and sigmin),
    output an integer z according to the discrete
    Gaussian distribution D_{Z, mu, sigma}.

    Params:
    - the center mu
    - the standard deviation sigma
    - a scaling factor sigmin(different between security levels of Falcon)
    - optional: the randomness source rng (default: OsRng)

    The inputs MUST verify 1 < sigmin < sigma < MAX_SIGMA.

    Returns:
    - a sample z from the distribution D_{Z, mu, sigma}.
*/

fn sampler_z<R: RngCore>(mu: f64, sigmin: f64, sigma: f64, rng: &mut R) -> i64 {
    let s = f64::floor(mu);
    let r = mu - s;
    let ccs = sigmin / sigma;
    let dss = 1.0 / (2.0 * sigma * sigma);

    loop {
        let z0 = base_sampler(rng);
        let b = rng.gen::<u8>() & 1;
        let z = b as i64 + (2 * b as i64 - 1) * z0 as i64;
        let z_minus_r = z as f64 - r;

        let mut x = (z_minus_r * z_minus_r) * dss;

        x -= (z0 * z0) as f64 * INV_2SIGMA2;

        if berexp(x, ccs, rng) {
            return z + s as i64;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hex_literal::hex;
    use std::io::ErrorKind;

    struct MockRng {
        bytes: Vec<u8>,
        index: usize,
    }

    impl MockRng {
        fn new(bytes: Vec<u8>) -> Self {
            MockRng {
                bytes: bytes,
                index: 0,
            }
        }

        fn next_u8(&mut self) -> u8 {
            let res = self.bytes[self.index];
            self.index += 1;

            res
        }
    }

    impl RngCore for MockRng {
        //Get next u32, pads with zeroes if self.bytes % 4 > 0
        fn next_u32(&mut self) -> u32 {
            if self.index >= self.bytes.len() {
                return 0; // Return 0 if there are no more bytes to read
            }

            let mut bytes = [0u8; 4];
            let len = std::cmp::min(4, self.bytes.len() - self.index);
            bytes[..len].copy_from_slice(&self.bytes[self.index..self.index + len]);
            self.index += len;

            u32::from_le_bytes(bytes)
        }

        fn next_u64(&mut self) -> u64 {
            let res =
                u64::from_le_bytes(self.bytes[self.index..self.index + 8].try_into().unwrap());
            self.index += 8;

            res
        }

        // Insecure implementation used for the purpose of mocking RNG
        fn fill_bytes(&mut self, dest: &mut [u8]) {
            let remaining = self.bytes.len().checked_sub(self.index).unwrap_or(0);
            if dest.len() > remaining {
                self.index = 0;
                panic!(
                    "MockRng fill_bytes() does not accept a dest longer than the remaining bytes"
                );
            }

            for i in dest.iter_mut() {
                *i = self.next_u8();
            }

            self.index = self
                .index
                .checked_add(dest.len())
                .unwrap_or(self.bytes.len());
        }

        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
            if dest.len() >= self.bytes.len() - self.index {
                eprintln!("MockRng fill_bytes() does not accept a dest longer than self.bytes");
                let error = std::io::Error::new(ErrorKind::Other, "Not enough bytes in MockRng");
                return Err(rand::Error::new(error));
            }

            for i in dest.iter_mut() {
                *i = self.next_u8();
            }
            self.index += dest.len();

            Ok(())
        }
    }

    #[derive(Debug)]
    struct Kat {
        mu: f64,
        sigma_prime: f64,
        random_bytes: &'static [u8],
        output_z: i32,
    }

    #[test]
    fn test_samplerz_kats() {
        /*
            Test vectors for SamplerZ as defined in the Falcon specification
            table 3.2 [https://falcon-sign.info/falcon.pdf]
        */
        let sigma_min: f64 = 1.277_833_697;
        let kats: [Kat; 16] = [
            Kat {
                sigma_prime: 1.7037990414754918,
                mu: -91.90471153063714,
                random_bytes: &hex!(
                    "0fc5442ff043d66e91d1ea
                    cac64ea5450a22941edc6c"
                ),
                output_z: -92,
            },
            Kat {
                sigma_prime: 1.7037990414754918,
                mu: -8.322564895434937,
                random_bytes: &hex!(
                    "f4da0f8d8444d1a77265c2
                    ef6f98bbbb4bee7db8d9b3"
                ),
                output_z: -8,
            },
            Kat {
                sigma_prime: 1.7035823083824078,
                mu: -19.096516109216804,
                random_bytes: &hex!(
                    "db47f6d7fb9b19f25c36d6
                    b9334d477a8bc0be68145d"
                ),
                output_z: -20,
            },
            Kat {
                sigma_prime: 1.7035823083824078,
                mu: -11.335543982423326,
                random_bytes: &hex!(
                    "ae41b4f5209665c74d00dc
                    c1a8168a7bb516b3190cb4
                    2c1ded26cd52aed770eca7
                    dd334e0547bcc3c163ce0b"
                ),
                output_z: -12,
            },
            Kat {
                sigma_prime: 1.6984647769450156,
                mu: 7.9386734193997555,
                random_bytes: &hex!(
                    "31054166c1012780c603ae
                    9b833cec73f2f41ca5807c
                    c89c92158834632f9b1555"
                ),
                output_z: 8,
            },
            Kat {
                sigma_prime: 1.6984647769450156,
                mu: -28.990850086867255,
                random_bytes: &hex!("737e9d68a50a06dbbc6477"),
                output_z: -30,
            },
            Kat {
                sigma_prime: 1.6980782114808988,
                mu: -9.071257914091655,
                random_bytes: &hex!("a98ddd14bf0bf22061d632"),
                output_z: -10,
            },
            Kat {
                sigma_prime: 1.6980782114808988,
                mu: -43.88754568839566,
                random_bytes: &hex!("3cbf6818a68f7ab9991514"),
                output_z: -41,
            },
            Kat {
                sigma_prime: 1.7010983419195522,
                mu: -58.17435547946095,
                random_bytes: &hex!(
                    "6f8633f5bfa5d26848668e
                    3d5ddd46958e97630410587c"
                ),
                output_z: -61,
            },
            Kat {
                sigma_prime: 1.7010983419195522,
                mu: -43.58664906684732,
                random_bytes: &hex!(
                    "272bc6c25f5c5ee53f83c4
                    3a361fbc7cc91dc783e20a"
                ),
                output_z: -46,
            },
            Kat {
                sigma_prime: 1.7009387219711465,
                mu: -34.70565203313315,
                random_bytes: &hex!(
                    "45443c59574c2c3b07e2e1
                    d9071e6d133dbe32754b0a"
                ),
                output_z: -34,
            },
            Kat {
                sigma_prime: 1.7009387219711465,
                mu: -44.36009577368896,
                random_bytes: &hex!(
                    "6ac116ed60c258e2cbaeab
                    728c4823e6da36e18d08da
                    5d0cc104e21cc7fd1f5ca8
                    d9dbb675266c928448059e"
                ),
                output_z: -44,
            },
            Kat {
                sigma_prime: 1.6958406126012802,
                mu: -21.783037079346236,
                random_bytes: &hex!("68163bc1e2cbf3e18e7426"),
                output_z: -23,
            },
            Kat {
                sigma_prime: 1.6958406126012802,
                mu: -39.68827784633828,
                random_bytes: &hex!("d6a1b51d76222a705a0259"),
                output_z: -40,
            },
            Kat {
                sigma_prime: 1.6955259305261838,
                mu: -18.488607061056847,
                random_bytes: &hex!(
                    "f0523bfaa8a394bf4ea5c1
                    0f842366fde286d6a30803"
                ),
                output_z: -22,
            },
            Kat {
                sigma_prime: 1.6955259305261838,
                mu: -48.39610939101591,
                random_bytes: &hex!(
                    "87bd87e63374cee62127fc
                    6931104aab64f136a0485b"
                ),
                output_z: -50,
            },
        ];

        let mut count = 1;
        for kat in kats {
            let bytes: &[u8] = kat.random_bytes;
            let mut rng = MockRng::new((bytes).to_vec());

            println!("Kat {count}");
            println!("bytes: {bytes:?}, bytes len: {}\n", bytes.len());

            let output = sampler_z(kat.mu, sigma_min, kat.sigma_prime, &mut rng);
            let expected_output = kat.output_z;
            assert_eq!(
                output, kat.output_z as i64,
                "Output was {output}, expected {expected_output}"
            );

            count += 1;
        }
    }
}
