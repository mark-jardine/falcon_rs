use num_bigint::{BigInt};

// Q is the modulus used in Falcon (12 * 1024 + 1)
pub const Q: u16 = 12 * 1024 + 1;

// This function raises base to the power of exponent, modulo modulus.
pub fn mod_pow(base: &BigInt, exponent: &BigInt, modulus: &BigInt) -> BigInt{
    BigInt::new(base.sign(), base.iter_u32_digits().collect()).modpow(exponent, modulus)
}

// Computes the modular inverse of number modulo modulus.
// It returns None if the inverse does not exist.
pub fn mod_inverse(number: &BigInt, modulus: &BigInt) -> Option<BigInt> {
    let (g, x, _) = extended_gcd(number, modulus);

    if g == BigInt::from(1) {
        // Ensure the inverse is positive
        Some((x % modulus + modulus) % modulus)
    } else {
        // Inverse does not exist
        None
    }
}

// Helper function which calculates the greatest common divisor
fn extended_gcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    if *a == BigInt::from(0) {
        (b.clone(), BigInt::from(0), BigInt::from(1))
    } else {
        let (g, x, y) = extended_gcd(&(b % a), a);
        (g, y - (b / a) * &x, x)
    }
}

// todo: determine if this is needed
// Finds the square root of number modulo modulus.
// pub fn mod_sqrt(number: &BigInt, modulus: &BigInt) -> BigInt{
//     number.sqrt().modpow(&BigInt::new(Sign::Plus, vec![1]), modulus)
// }

mod tests{
    use super::*;
    use num_bigint::ToBigInt;

    #[test]
    fn test_mod_pow() {
        assert_eq!(
            mod_pow(&2.to_bigint().unwrap(), &3.to_bigint().unwrap(), &5.to_bigint().unwrap()),
            3.to_bigint().unwrap()
        );
        // Add more test cases...
    }

    #[test]
    fn test_mod_inverse() {
        assert_eq!(
            mod_inverse(&3.to_bigint().unwrap(), &11.to_bigint().unwrap()),
            Some(4.to_bigint().unwrap())
        );
        assert_eq!(
            mod_inverse(&2.to_bigint().unwrap(), &4.to_bigint().unwrap()),
            None
        );
        // Add more test cases...
    }

    // Determine desired behaviour
    // #[test]
    // fn test_mod_sqrt() {
    //
    // }

}