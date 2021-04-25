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
    pub fn new(x: K, y: K) -> Result<Self, Error<K>> {
        if !x.is_positive() || !y.is_positive() {
            Err(Error::InvalidChromacity(x, y))
        } else {
            Ok(Self(x, y))
        }
    }

    pub unsafe fn new_unchecked(x: K, y: K) -> Self { Self(x, y) }
}

impl<K: matrix::Scalar> Chromacity<K>
where
    for<'x> &'x K: num_traits::RefNum<K>,
{
    /// Returns XYZ coordinates of a colour with given chromacity and luminosity
    /// (the Y coordinate) equal one.
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


pub use matrix::calculate as calculate_matrix;


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
            super::calculate_matrix(&white(f).to_xyz(), &primaries(f))
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
            super::calculate_matrix(&white(f).to_xyz(), &primaries(f))
        );
    }

    #[test]
    fn test_calculate_ratio() { run_ratio_test(&new_ratio); }

    #[test]
    fn test_calculate_big_ratio() { run_ratio_test(&new_big_ratio); }
}
