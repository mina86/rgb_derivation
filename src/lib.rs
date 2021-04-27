//! RGB colour system derivation routines.
//!
//! Functions for deriving RGB→XYZ and XYZ→RGB conversion matrices for given RGB
//! colour system (such as sRGB colour space).  The calculations are performed
//! from the definition of such system provided in the form of chromacities of
//! the reference white point and the red, green and blue primaries.
//! Alternatively, constructions from XYZ coordinates of primaries is also
//! available.
//!
//! The crate supports arithmetic with rational and big integer types such that
//! the calculations can be performed without any loss of precision if desired.
//! So long as a type implements the four basic arithmetic operations, it can be
//! used with this library.  For example, `f32`, `num::Rational64` and
//! `num::BigRational` can all be used.
//!
//! # Example
//!
//! ```
//! type Scalar = num::BigRational;
//! type Chromacity = rgb_derivation::Chromacity<Scalar>;
//!
//! fn scalar(numer: i64, denom: i64) -> num::BigRational {
//!     num::BigRational::new(numer.into(), denom.into())
//! }
//!
//! fn chromacity(x: (i64, i64), y: (i64, i64)) -> Chromacity {
//!     Chromacity::new(scalar(x.0, x.1), scalar(y.0, y.1)).unwrap()
//! }
//!
//! let white_xy = chromacity((31, 100), (33, 100));
//! let primaries_xy = [
//!     chromacity((64, 100), (33, 100)),
//!     chromacity((30, 100), (60, 100)),
//!     chromacity((15, 100), (6, 100)),
//! ];
//!
//! let white_xyz = white_xy.to_xyz();
//! let matrix =
//!     rgb_derivation::matrix::calculate(&white_xyz, &primaries_xy).unwrap();
//! let inverse = rgb_derivation::matrix::inversed_copy(&matrix).unwrap();
//! let primaries_xyz = rgb_derivation::matrix::transposed_copy(&matrix);
//!
//! assert_eq!([scalar(31, 33), scalar(1, 1), scalar(12, 11)], white_xyz);
//! assert_eq!([
//!     [scalar(1088, 2739), scalar(17, 83), scalar(17, 913)],
//!     [scalar(  30,   83), scalar(60, 83), scalar(10,  83)],
//!     [scalar(  15,   83), scalar( 6, 83), scalar(79,  83)]
//! ], primaries_xyz);
//! assert_eq!([
//!     [scalar(1088, 2739), scalar(30, 83), scalar(15, 83)],
//!     [scalar(  17,   83), scalar(60, 83), scalar( 6, 83)],
//!     [scalar(  17,  913), scalar(10, 83), scalar(79, 83)]
//! ], matrix);
//! assert_eq!([
//!     [scalar( 286,  85), scalar(-407,  255), scalar(-44,  85)],
//!     [scalar(-863, 900), scalar(5011, 2700), scalar( 37, 900)],
//!     [scalar(   1,  18), scalar( -11,   54), scalar( 19,  18)]
//! ], inverse);
//! ```
//!
//! (Note: if you need matrices for the sRGB colour space, the [`srgb`
//! crate](https://crates.io/crates/srgb) provides them along with gamma
//! functions needed to properly handle sRGB)

pub mod matrix;


/// Possible errors which can occur when performing calculations.
#[derive(PartialEq, Eq, Debug)]
pub enum Error<K> {
    /// Error returned when trying to create [`Chromacity`] with either of the
    /// coordinates being a non-positive number.  The two arguments are the
    /// coordinates which caused the issue.
    InvalidChromacity(K, K),
    /// Error returned if provided reference white point has non-positive
    /// luminosity (the `Y` component).  The argument are the XYZ coordinates of
    /// the white point which caused the issue.
    InvalidWhitePoint([K; 3]),
    /// Error returned when XYZ coordinates of primaries are linearly dependent.
    /// That is, when one of the primaries is a linear combination of the other
    /// two.
    DegenerateMatrix,
}


/// A colour chromacity represented as `(x, y)` coordinates.
#[derive(Clone, Copy, Debug, PartialOrd, Ord, PartialEq, Eq)]
pub struct Chromacity<K>(K, K);

impl<K> Chromacity<K> {
    pub fn x(&self) -> &K { &self.0 }
    pub fn y(&self) -> &K { &self.1 }
}

impl<K: num_traits::Signed> Chromacity<K> {
    /// Constructs new Chromacity from given (x, y) coordinates.
    ///
    /// Returns an error if either of the coordinate is non-positive.
    pub fn new(x: K, y: K) -> Result<Self, Error<K>> {
        if !x.is_positive() || !y.is_positive() {
            Err(Error::InvalidChromacity(x, y))
        } else {
            Ok(Self(x, y))
        }
    }

    /// Constructs new Chromacity from given (x, y) coordinates.
    ///
    /// Does not check whether the coordinates are positive.  If they aren’t,
    /// other methods (e.g. [`Chromacity::to_xyz`] may result in undefined
    /// behaviour.
    pub unsafe fn new_unchecked(x: K, y: K) -> Self { Self(x, y) }
}

impl<K: matrix::Scalar + num_traits::FromPrimitive> Chromacity<K>
where
    for<'x> &'x K: num_traits::RefNum<K> {
    /// Returns (u′, v′) coordinates of the chromacity.
    ///
    /// **Panics** if `K` cannot be constructed from small primitive integers,
    /// i.e. when `K::from_i32` returns `None`.
    ///
    /// # Example
    ///
    /// ```
    /// use rgb_derivation::*;
    ///
    /// let one_third = num::rational::Ratio::new(1i64, 3i64);
    /// let uv = Chromacity::new(one_third, one_third).unwrap().to_uv();
    /// assert_eq!((
    ///     num::rational::Ratio::new(4, 19),
    ///     num::rational::Ratio::new(9, 19)
    /// ), uv);
    /// ```
    pub fn to_uv(&self) -> (K, K) {
        let int = |n| K::from_i32(n).unwrap();
        let denom = int(-2) * self.x() + int(12) * self.y() + int(3);
        let u = int(4) * self.x() / &denom;
        let v = int(9) * self.y() / denom;
        (u, v)
    }
}

impl<K: matrix::Scalar> Chromacity<K>
where
    for<'x> &'x K: num_traits::RefNum<K>,
{
    /// Returns XYZ coordinates of a colour with given chromacity.  Assumes
    /// luminosity (the Y coordinate) equal one.
    ///
    /// # Example
    ///
    /// ```
    /// use rgb_derivation::*;
    ///
    /// let one = num::rational::Ratio::new(1i64, 1i64);
    /// let one_third = num::rational::Ratio::new(1i64, 3i64);
    /// let e = Chromacity::new(one_third, one_third).unwrap().to_xyz();
    /// assert_eq!([one, one, one], e);
    /// ```
    pub fn to_xyz(&self) -> [K; 3] {
        let (x, y) = (self.x(), self.y());
        let uc_x = x / y;
        let uc_z = (K::one() - x - y) / y;
        [uc_x, K::one(), uc_z]
    }
}


#[cfg(test)]
pub(crate) mod test {
    type Ratio = (i64, i64);

    pub(crate) fn new_float(num: Ratio) -> f32 {
        (num.0 as f64 / num.1 as f64) as f32
    }

    pub(crate) fn new_ratio(num: Ratio) -> num::rational::Ratio<i128> {
        num::rational::Ratio::new(num.0 as i128, num.1 as i128)
    }

    pub(crate) fn new_big_ratio(num: Ratio) -> num::BigRational {
        num::rational::Ratio::new(num.0.into(), num.1.into())
    }

    pub(crate) fn chromacity<K>(
        f: &impl Fn(Ratio) -> K,
        x: Ratio,
        y: Ratio,
    ) -> super::Chromacity<K>
    where
        K: std::fmt::Debug + num_traits::Signed, {
        super::Chromacity::new(f(x), f(y)).unwrap()
    }


    fn white<K>(f: &impl Fn((i64, i64)) -> K) -> super::Chromacity<K>
    where
        K: std::fmt::Debug + num_traits::Signed, {
        chromacity(f, (312713, 1000000), (41127, 125000))
    }

    fn primaries<K>(f: &impl Fn((i64, i64)) -> K) -> [super::Chromacity<K>; 3]
    where
        K: std::fmt::Debug + num_traits::Signed, {
        [
            chromacity(f, (64, 100), (33, 100)),
            chromacity(f, (30, 100), (60, 100)),
            chromacity(f, (15, 100), (06, 100)),
        ]
    }

    #[test]
    fn test_xyz() {
        let f = &new_float;
        assert_eq!([0.9504492, 1.0, 1.0889165], white(f).to_xyz());

        let f = &new_ratio;
        assert_eq!(
            [f((312713, 329016)), f((1, 1)), f((358271, 329016))],
            white(f).to_xyz()
        );

        let f = &new_big_ratio;
        assert_eq!(
            [f((312713, 329016)), f((1, 1)), f((358271, 329016))],
            white(f).to_xyz()
        );
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
            super::matrix::calculate(&white(f).to_xyz(), &primaries(f))
        );
    }

    fn run_ratio_test<K>(f: &impl Fn((i64, i64)) -> K)
    where
        K: super::matrix::Scalar + num_traits::Signed + std::fmt::Debug,
        for<'x> &'x K: num_traits::RefNum<K>, {
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
            super::matrix::calculate(&white(f).to_xyz(), &primaries(f))
        );
    }

    #[test]
    fn test_calculate_ratio() { run_ratio_test(&new_ratio); }

    #[test]
    fn test_calculate_big_ratio() { run_ratio_test(&new_big_ratio); }
}
