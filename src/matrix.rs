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

//! Functions for calculating RGB↔XYZ conversion matrices and performing basic
//! matrix manipulation.
//!
//! Specifically, [`calculate`] generates the RGB→XYZ change of basis matrix
//! from chromacicities of the reference white point and red, green and blue
//! primary colours.  Inversing that matrix with [`inversed_copy`] constructs
//! change of basis in the opposite direction and transposing with
//! [`transposed_copy`] results in a matrix whose rows are XYZ coordinates of
//! the primary colours.


/// Trait for scalar type used in calculations.
///
/// An implementation of this trait is provided to all types which satisfy this
/// traits bounds.
pub trait Scalar:
    Clone + num_traits::NumRef + num_traits::NumAssignRef + num_traits::Signed {
}

impl<T> Scalar for T where
    T: Clone
        + core::ops::Neg<Output = Self>
        + num_traits::NumRef
        + num_traits::NumAssignRef
        + num_traits::Signed
{
}


/// A two-dimensional 3×3 array.
///
/// The crate assumes a row-major order for the matrix, i.e. the first index
/// specifies the row and the second index specifies column within that row.
pub type Matrix<K> = [[K; 3]; 3];


/// Calculates change of basis matrix for moving from linear RGB to XYZ colour
/// spaces.
///
/// The matrix is calculated from XYZ coordinates of a reference white point and
/// chromacicities of the three primary colours (red, green and blue).  (Note that
/// [`crate::Chromaticity::to_xyz`] function allows conversion from
/// chromaticity to XYZ coordinates thus the function may be used when only x and
/// y coordinates of the white point are known).
///
/// The result is a three-by-three matrix M such that multiplying it by
/// a column-vector representing a colour in linear RGB space results in
/// a column-vector representing the same colour in XYZ coordinates.
///
/// To get the change of basis matrix for moving in the other direction
/// (i.e. from XYZ colour space to linear RGB space) simply inverse the result.
/// This scan be done with [`inversed_copy`] function.
///
/// Finally, the columns of the result are XYZ coordinates of the primaries.  To
/// get the primaries oriented in rows [`transposed_copy`] function can be used.
///
/// # Example
///
/// ```
/// use rgb_derivation::{Chromaticity, matrix};
///
/// let white = Chromaticity::new(1.0_f32 / 3.0, 1.0 / 3.0).unwrap();
/// let primaries = [
///     Chromaticity::new(0.735_f32, 0.265_f32).unwrap(),
///     Chromaticity::new(0.274_f32, 0.717_f32).unwrap(),
///     Chromaticity::new(0.167_f32, 0.009_f32).unwrap(),
/// ];
///
/// let matrix    = matrix::calculate(&white.to_xyz(), &primaries).unwrap();
/// let inverse   = matrix::inversed_copy(&matrix).unwrap();
/// let primaries = matrix::transposed_copy(&matrix);
///
/// assert_eq!([
///     [0.4887181,  0.31068033,  0.20060167],
///     [0.17620447, 0.8129847,   0.010810869],
///     [0.0,        0.010204833, 0.98979515],
/// ], matrix);
/// assert_eq!([
///     [ 2.3706737,   -0.9000402,  -0.47063363],
///     [-0.5138849,    1.4253035,   0.08858136],
///     [ 0.005298177, -0.014694944, 1.0093968],
/// ], inverse);
/// assert_eq!([
///     [0.4887181,  0.17620447,  0.0],
///     [0.31068033, 0.8129847,   0.010204833],
///     [0.20060167, 0.010810869, 0.98979515],
/// ], primaries);
/// ```
pub fn calculate<K: Scalar>(
    white: &[K; 3],
    primaries: &[crate::Chromaticity<K>; 3],
) -> Result<Matrix<K>, crate::Error<K>>
where
    for<'x> &'x K: num_traits::RefNum<K>, {
    if !white[1].is_positive() {
        return Err(crate::Error::InvalidWhitePoint(white.clone()));
    }

    // Calculate the transformation matrix as per
    // https://mina86.com/2019/srgb-xyz-matrix/

    // M' = [[R_X/R_Y 1 R_Z/R_Y] [G_X/G_Y 1 G_Z/G_Y] [B_X/B_Y 1 B_Z/B_Y]]^T
    let mut mp = transposed(make_vector(|i| primaries[i].to_xyz()));

    // Y = M′⁻¹ ✕ W
    let inverse_m_prime = inversed_copy(&mp)?;
    let y_fn = |i| dot_product(&inverse_m_prime[i], &white);

    // M = M′ ✕ diag(Y)
    for col in 0..3 {
        let y = y_fn(col);
        for row in 0..3 {
            mp[row][col] *= &y;
        }
    }

    Ok(mp)
}


/// Transposes a 3×3 matrix.  Consumes the argument and returns a new matrix.
///
/// # Example
///
/// ```
/// let primaries = [
///     [0.431, 0.222, 0.020],
///     [0.342, 0.707, 0.130],
///     [0.178, 0.071, 0.939]
/// ];
/// assert_eq!([
///    [0.431, 0.342, 0.178],
///    [0.222, 0.707, 0.071],
///    [0.020, 0.130, 0.939],
/// ], rgb_derivation::matrix::transposed(primaries));
/// ```
pub fn transposed<K>(mut matrix: Matrix<K>) -> Matrix<K> {
    let m = matrix.as_mut_ptr();
    for (i, j) in [(0, 1), (0, 2), (1, 2)].iter().copied() {
        // SAFETY: In each step, i != j therefore matrix[i] and matrix[j]
        // are different rows.
        let (a, b) =
            unsafe { (&mut *m.offset(i as isize), &mut *m.offset(j as isize)) };
        core::mem::swap(&mut a[j], &mut b[i]);
    }
    matrix
}

/// Transposes a 3×3 matrix.  Constructs a new matrix and returns it.
///
/// This is equivalent to [`transposed`] except that it doesn’t consume the
/// argument and returns a new object.
pub fn transposed_copy<K: Clone>(matrix: &Matrix<K>) -> Matrix<K> {
    make_matrix(|i, j| matrix[j][i].clone())
}


/// Returns inversion of a 3✕3 matrix M, i.e. M⁻¹.
///
/// Returns an error if the matrix is non-invertible, i.e. if it’s determinant
/// is zero.
///
/// # Example
///
/// ```
/// let matrix: [[f32; 3]; 3] = [
///     [0.488, 0.310, 0.200],
///     [0.176, 0.812, 0.010],
///     [0.000, 0.010, 0.989],
/// ];
/// assert_eq!([
///     [ 2.3739555,    -0.90051305, -0.4709666],
///     [-0.5146161,     1.4268899,   0.08964035],
///     [ 0.0052033975, -0.014427602, 1.010216],
/// ], rgb_derivation::matrix::inversed_copy(&matrix).unwrap());
/// ```
pub fn inversed_copy<K>(
    matrix: &Matrix<K>,
) -> Result<Matrix<K>, crate::Error<K>>
where
    K: Scalar,
    for<'x> &'x K: num_traits::RefNum<K>, {
    let mut comatrix_transposed =
        make_matrix(|row, col| cofactor(matrix, col, row));

    // https://en.wikipedia.org/wiki/Minor_(linear_algebra)#Cofactor_expansion_of_the_determinant
    // Because we transposed the comatrix when we created it, we need to
    // calculate dot product of the first row of the matrix and first *column*
    // of the comatrix.
    let det: K = dot_product_with_column(&matrix[0], &comatrix_transposed, 0);
    if det.is_zero() {
        return Err(crate::Error::DegenerateMatrix);
    }

    // https://en.wikipedia.org/wiki/Minor_(linear_algebra)#Inverse_of_a_matrix
    // We’ve already transposed comatrix so now all we have to do is just divide
    // it by the Scalar.
    for row in 0..3 {
        for col in 0..3 {
            comatrix_transposed[row][col] /= &det;
        }
    }
    Ok(comatrix_transposed)
}



/// Constructs a new 3-element array by applying a function to its indices.
fn make_vector<T>(f: impl Fn(usize) -> T) -> [T; 3] { [f(0), f(1), f(2)] }

/// Constructs a new 2-dimensional 3×3 array by applying a function to all of
/// its indices.
fn make_matrix<T>(f: impl Fn(usize, usize) -> T) -> [[T; 3]; 3] {
    make_vector(|r| make_vector(|c| f(r, c)))
}

/// Calculates a dot product of two 3-element vectors.
fn dot_product<K: Scalar>(a: &[K; 3], b: &[K; 3]) -> K
where
    for<'x> &'x K: num_traits::RefNum<K>, {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

/// Calculates a dot product of a 3-element vectors and a column in a matrix.
fn dot_product_with_column<K>(a: &[K; 3], b: &Matrix<K>, col: usize) -> K
where
    K: Scalar,
    for<'x> &'x K: num_traits::RefNum<K>, {
    &a[0] * &b[0][col] + &a[1] * &b[1][col] + &a[2] * &b[2][col]
}


/// Returns a cofactor of a 3✕3 matrix M, i.e. C_{row,col}.
fn cofactor<K: Scalar>(matrix: &Matrix<K>, row: usize, col: usize) -> K
where
    for<'x> &'x K: num_traits::RefNum<K>, {
    let rr = ((row == 0) as usize, 2 - (row == 2) as usize);
    let cc = ((col == 0) as usize, 2 - (col == 2) as usize);
    let ad = &matrix[rr.0][cc.0] * &matrix[rr.1][cc.1];
    let bc = &matrix[rr.1][cc.0] * &matrix[rr.0][cc.1];
    let minor = ad - bc;
    if (row ^ col) & 1 == 0 {
        minor
    } else {
        -minor
    }
}


#[test]
fn test_transpose() {
    let matrix: Matrix<u64> = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
    assert_eq!([[1, 4, 7], [2, 5, 8], [3, 6, 9]], transposed_copy(&matrix));

    let matrix: Matrix<u64> = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
    assert_eq!([[1, 4, 7], [2, 5, 8], [3, 6, 9]], transposed(matrix));
}

#[test]
fn test_inverse_floats() {
    assert_eq!(
        Ok([
            [3.240812809834622, -1.5373086942720335, -0.49858660478241557],
            [-0.9692430382347864, 1.8759663312198533, 0.04155504934405438],
            [
                0.05563834281000593,
                -0.20400734898293651,
                1.0571294977107015
            ]
        ]),
        inversed_copy(&[
            [0.4124108, 0.35758457, 0.18045382],
            [0.21264932, 0.71516913, 0.07218152],
            [0.019331757, 0.119194806, 0.9503901],
        ])
    );
}

#[cfg(test)]
fn run_inverse_ratio_test<K>(f: &impl Fn((i64, i64)) -> K)
where
    K: Scalar + core::fmt::Debug,
    for<'x> &'x K: num_traits::RefNum<K>, {
    assert_eq!(
        Ok([
            [
                f((4277208, 1319795)),
                f((-2028932, 1319795)),
                f((-658032, 1319795))
            ],
            [
                f((-70985202, 73237775)),
                f((137391598, 73237775)),
                f((3043398, 73237775))
            ],
            [
                f((164508, 2956735)),
                f((-603196, 2956735)),
                f((3125652, 2956735))
            ]
        ]),
        inversed_copy(&[
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
        ])
    );
}

#[test]
fn test_inverses_ratio() { run_inverse_ratio_test(&crate::test::new_ratio); }

#[test]
fn test_inverses_big_ratio() {
    run_inverse_ratio_test(&crate::test::new_big_ratio);
}
