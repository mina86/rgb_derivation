/* This file is part of rgb_derivation crate.
 * Copyright 2021 by Michał Nazarewicz <mina86@mina86.com>
 *
 * rgb_derivation crate is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at your
 * option) any later version.
 *
 * rgb_derivation crate is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * rgb_derivation crate.  If not, see <http://www.gnu.org/licenses/>. */

#![no_std]
#![doc = include_str!("../README.md")]

pub mod matrix;


/// Possible errors which can occur when performing calculations.
#[derive(PartialEq, Eq, Debug)]
pub enum Error<K> {
    /// Error returned when trying to create [`Chromaticity`] with either of the
    /// coordinates being a non-positive number.  The two arguments are the
    /// coordinates which caused the issue.
    InvalidChromaticity(K, K),
    /// Error returned if provided reference white point has non-positive
    /// luminosity (the `Y` component).  The argument are the XYZ coordinates of
    /// the white point which caused the issue.
    InvalidWhitePoint([K; 3]),
    /// Error returned when XYZ coordinates of primaries are linearly dependent.
    /// That is, when one of the primaries is a linear combination of the other
    /// two.
    DegenerateMatrix,
}


/// A colour chromaticity represented as `(x, y)` coordinates.
#[derive(Clone, Copy, Debug, PartialOrd, Ord, PartialEq, Eq)]
pub struct Chromaticity<K>(K, K);

impl<K> Chromaticity<K> {
    pub fn x(&self) -> &K { &self.0 }
    pub fn y(&self) -> &K { &self.1 }

    /// Deconstructs the object into `(x, y)` pair.
    pub fn into_xy(self) -> (K, K) { (self.0, self.1) }
}

impl<K: num_traits::Signed> Chromaticity<K> {
    /// Constructs new Chromaticity from given (x, y) coordinates.
    ///
    /// Returns an error if either of the coordinate is non-positive.
    pub fn new(x: K, y: K) -> Result<Self, Error<K>> {
        if !x.is_positive() || !y.is_positive() {
            Err(Error::InvalidChromaticity(x, y))
        } else {
            Ok(Self(x, y))
        }
    }

    /// Constructs new Chromaticity from given (x, y) coordinates.
    ///
    /// # Safety
    ///
    /// Does not check whether the coordinates are positive.  If they aren’t,
    /// other methods (e.g. [`Chromaticity::into_xyz`] may result in undefined
    /// behaviour.
    pub unsafe fn new_unchecked(x: K, y: K) -> Self { Self(x, y) }
}

impl<K: matrix::Scalar> Chromaticity<K>
where
    for<'x> &'x K: num_traits::NumOps<&'x K, K>,
{
    /// Returns XYZ coordinates of a colour with given chromaticity.  Assumes
    /// luminosity (the Y coordinate) equal one.
    ///
    /// # Example
    ///
    /// ```
    /// use rgb_derivation::Chromaticity;
    ///
    /// let one = num::rational::Ratio::new(1i64, 1i64);
    /// let one_third = num::rational::Ratio::new(1i64, 3i64);
    /// let e = Chromaticity::new(one_third, one_third).unwrap().into_xyz();
    /// assert_eq!([one, one, one], e);
    /// ```
    pub fn to_xyz(&self) -> [K; 3] {
        let (x, y) = (self.x(), self.y());
        let uc_x = x / y;
        let uc_z = (K::one() - x - y) / y;
        [uc_x, K::one(), uc_z]
    }
}

impl<K: matrix::Scalar> Chromaticity<K> {
    /// Returns XYZ coordinates of a colour with given chromaticity.  Assumes
    /// luminosity (the Y coordinate) equal one.
    ///
    /// # Example
    ///
    /// ```
    /// use rgb_derivation::Chromaticity;
    ///
    /// let one = num::rational::Ratio::new(1i64, 1i64);
    /// let one_third = num::rational::Ratio::new(1i64, 3i64);
    /// let e = Chromaticity::new(one_third, one_third).unwrap().into_xyz();
    /// assert_eq!([one, one, one], e);
    /// ```
    pub fn into_xyz(self) -> [K; 3] {
        let (x, y) = self.into_xy();
        let uc_z = (K::one() - &x - &y) / &y;
        let uc_x = x / y;
        [uc_x, K::one(), uc_z]
    }
}


#[cfg(test)]
pub(crate) mod test {
    type Ratio = (i32, i32);
    pub type Rational = num::rational::Ratio<i128>;

    pub fn float(num: Ratio) -> f32 { (num.0 as f64 / num.1 as f64) as f32 }

    pub const fn ratio(num: Ratio) -> Rational {
        assert!(num.1 != 0);
        num::rational::Ratio::new_raw(num.0 as i128, num.1 as i128)
    }

    pub fn chromaticity<K: core::fmt::Debug + num_traits::Signed>(
        f: &impl Fn(Ratio) -> K,
        x: Ratio,
        y: Ratio,
    ) -> super::Chromaticity<K> {
        super::Chromaticity::new(f(x), f(y)).unwrap()
    }


    pub fn white<K: core::fmt::Debug + num_traits::Signed>(
        f: &impl Fn((i32, i32)) -> K,
    ) -> super::Chromaticity<K> {
        chromaticity(f, (312713, 1000000), (41127, 125000))
    }

    pub fn primaries<K: core::fmt::Debug + num_traits::Signed>(
        f: &impl Fn((i32, i32)) -> K,
    ) -> [super::Chromaticity<K>; 3] {
        [
            chromaticity(f, (64, 100), (33, 100)),
            chromaticity(f, (30, 100), (60, 100)),
            chromaticity(f, (15, 100), (06, 100)),
        ]
    }

    pub const M_F32: [[f32; 3]; 3] = [
        [0.4124108, 0.3575846, 0.18045382],
        [0.21264933, 0.7151692, 0.07218152],
        [0.019331757, 0.11919481, 0.9503901],
    ];

    pub const M_Q: [[Rational; 3]; 3] = [
        [
            ratio((4223344, 10240623)),
            ratio((14647555, 40962492)),
            ratio((14783675, 81924984)),
        ],
        [
            ratio((2903549, 13654164)),
            ratio((14647555, 20481246)),
            ratio((2956735, 40962492)),
        ],
        [
            ratio((263959, 13654164)),
            ratio((14647555, 122887476)),
            ratio((233582065, 245774952)),
        ],
    ];
    pub const M_INV_Q: [[Rational; 3]; 3] = [
        [
            ratio((4277208, 1319795)),
            ratio((-2028932, 1319795)),
            ratio((-658032, 1319795)),
        ],
        [
            ratio((-70985202, 73237775)),
            ratio((137391598, 73237775)),
            ratio((3043398, 73237775)),
        ],
        [
            ratio((164508, 2956735)),
            ratio((-603196, 2956735)),
            ratio((3125652, 2956735)),
        ],
    ];


    #[test]
    fn test_xyz() {
        let f = &float;
        let w = white(f);
        assert_eq!([0.9504492, 1.0, 1.0889165], w.to_xyz());
        assert_eq!([0.9504492, 1.0, 1.0889165], w.into_xyz());

        let f = &ratio;
        let w = white(f);
        let want = [f((312713, 329016)), f((1, 1)), f((358271, 329016))];
        assert_eq!(want, w.to_xyz());
        assert_eq!(want, w.into_xyz());
    }

    #[test]
    fn test_calculate_matrix_floats() {
        let f = &float;
        assert_eq!(
            Ok(M_F32),
            super::matrix::calculate(&white(f).into_xyz(), &primaries(f))
        );
    }

    #[test]
    fn test_calculate_ratio() {
        let f = &ratio;
        assert_eq!(
            Ok(M_Q),
            super::matrix::calculate(&white(f).into_xyz(), &primaries(f))
        );
    }
}
