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

    pub(crate) fn new_float(num: Ratio) -> f32 {
        (num.0 as f64 / num.1 as f64) as f32
    }

    pub(crate) fn new_ratio(num: Ratio) -> num::rational::Ratio<i128> {
        num::rational::Ratio::new(num.0 as i128, num.1 as i128)
    }

    fn chromaticity<K: core::fmt::Debug + num_traits::Signed>(
        f: &impl Fn(Ratio) -> K,
        x: Ratio,
        y: Ratio,
    ) -> super::Chromaticity<K> {
        super::Chromaticity::new(f(x), f(y)).unwrap()
    }


    fn white<K: core::fmt::Debug + num_traits::Signed>(
        f: &impl Fn((i32, i32)) -> K,
    ) -> super::Chromaticity<K> {
        chromaticity(f, (312713, 1000000), (41127, 125000))
    }

    fn primaries<K: core::fmt::Debug + num_traits::Signed>(
        f: &impl Fn((i32, i32)) -> K,
    ) -> [super::Chromaticity<K>; 3] {
        [
            chromaticity(f, (64, 100), (33, 100)),
            chromaticity(f, (30, 100), (60, 100)),
            chromaticity(f, (15, 100), (06, 100)),
        ]
    }

    #[test]
    fn test_xyz() {
        let f = &new_float;
        let w = white(f);
        assert_eq!([0.9504492, 1.0, 1.0889165], w.to_xyz());
        assert_eq!([0.9504492, 1.0, 1.0889165], w.into_xyz());

        let f = &new_ratio;
        let w = white(f);
        let want = [f((312713, 329016)), f((1, 1)), f((358271, 329016))];
        assert_eq!(want, w.to_xyz());
        assert_eq!(want, w.into_xyz());
    }

    #[test]
    fn test_calculate_matrix_floats() {
        let f = &new_float;
        assert_eq!(
            Ok([
                [0.4124108, 0.35758457, 0.18045382],
                [0.21264932, 0.71516913, 0.07218152],
                [0.019331757, 0.119194806, 0.9503901]
            ]),
            super::matrix::calculate(&white(f).into_xyz(), &primaries(f))
        );
    }

    #[test]
    fn test_calculate_ratio() {
        let f = &new_ratio;
        assert_eq!(
            Ok([
                [
                    f((4223344, 10240623)),
                    f((14647555, 40962492)),
                    f((14783675, 81924984))
                ],
                [
                    f((2903549, 13654164)),
                    f((14647555, 20481246)),
                    f((2956735, 40962492))
                ],
                [
                    f((263959, 13654164)),
                    f((14647555, 122887476)),
                    f((233582065, 245774952))
                ],
            ]),
            super::matrix::calculate(&white(f).into_xyz(), &primaries(f))
        );
    }
}
