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
    Clone + num_traits::NumRef + num_traits::NumAssignRef + num_traits::Signed
{
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
/// chromacicities of the three primary colours (red, green and blue).  (Note
/// that [`crate::Chromaticity::into_xyz`] function allows conversion from
/// chromaticity to XYZ coordinates thus the function may be used when only
/// x and y coordinates of the white point are known).
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
/// // CIE RGB
/// let white = Chromaticity::new(1.0_f32 / 3.0, 1.0 / 3.0).unwrap();
/// let primaries = [
///     Chromaticity::new(0.73474284_f32, 0.26525716_f32).unwrap(),
///     Chromaticity::new(0.27377903_f32, 0.7174777_f32).unwrap(),
///     Chromaticity::new(0.16655563_f32, 0.00891073_f32).unwrap(),
/// ];
///
/// let matrix  = matrix::calculate(&white.into_xyz(), &primaries).unwrap();
/// let inverse = matrix::inversed_copy(&matrix).unwrap();
///
/// assert_eq!([
///     [0.4900001,    0.31,        0.19999997],
///     [0.17690003,   0.8123999,   0.010700002],
///     [1.9875172e-8, 0.009900022, 0.99009985],
/// ], matrix);
/// assert_eq!([
///     [2.3644896,    -0.8965525,  -0.4679374],
///     [-0.51493526,   1.426333,    0.08860245],
///     [0.0051487973, -0.014261905, 1.0091132],
/// ], inverse);
/// ```
pub fn calculate<K: Scalar>(
    white: &[K; 3],
    primaries: &[crate::Chromaticity<K>; 3],
) -> Result<Matrix<K>, crate::Error<K>>
where
    for<'x> &'x K: num_traits::NumOps<&'x K, K>,
{
    if !white[1].is_positive() {
        return Err(crate::Error::InvalidWhitePoint(white.clone()));
    }

    // Calculate the transformation matrix as per
    // https://mina86.com/2019/srgb-xyz-matrix/

    // M' = [[R_X/R_Y 1 R_Z/R_Y] [G_X/G_Y 1 G_Z/G_Y] [B_X/B_Y 1 B_Z/B_Y]]^T
    let m_prime = transposed(make_vector(|i| primaries[i].to_xyz()));

    // Y = M'⁻¹ ✕ W = C'^T ✕ W / det(M').  By pulling det(M') out, we do 3
    // divisions rather than 9.  It also improves stability of the calculation
    // (as tested by calculating square error between M×M⁻¹ and identity
    // matrix).
    let (c_prime_transp, det) = inversed_internal(&m_prime)?;
    let y_fn = |i| dot_product(&c_prime_transp[i], white) / &det;

    // M = M' ✕ diag(Y)
    let mut m = m_prime;
    for col in 0..3 {
        let y = y_fn(col);
        for row in m.iter_mut() {
            row[col] *= &y;
        }
    }

    Ok(m)
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
    for (i, j) in [(0, 1), (0, 2), (1, 2)].iter().copied() {
        // Note: i < j
        let (a, b) = matrix.split_at_mut(j);
        core::mem::swap(&mut a[i][j], &mut b[0][i]);
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

#[test]
fn test_transpose() {
    let matrix: Matrix<u64> = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
    assert_eq!([[1, 4, 7], [2, 5, 8], [3, 6, 9]], transposed_copy(&matrix));
    assert_eq!([[1, 4, 7], [2, 5, 8], [3, 6, 9]], transposed(matrix));
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
pub fn inversed_copy<K: Scalar>(
    matrix: &Matrix<K>,
) -> Result<Matrix<K>, crate::Error<K>>
where
    for<'x> &'x K: num_traits::NumOps<&'x K, K>,
{
    inversed_internal(matrix).map(|(mut m, det)| {
        // https://en.wikipedia.org/wiki/Minor_(linear_algebra)#Inverse_of_a_matrix
        // We’ve already transposed comatrix so now all we have to do is divide
        // it by the determinant.
        for row in m.iter_mut() {
            for cell in row.iter_mut() {
                *cell /= &det;
            }
        }
        m
    })
}

/// For a given 3×3 matrix, returns `(com^T, 1/det)` pair where `com` is
/// comatrix of the matrix and `det` is its determinant.  To get inverse of the
/// matrix, multiply each element of `com^T` by `1/det`.
fn inversed_internal<K: Scalar>(
    matrix: &Matrix<K>,
) -> Result<(Matrix<K>, K), crate::Error<K>>
where
    for<'x> &'x K: num_traits::NumOps<&'x K, K>,
{
    let comatrix_transposed =
        make_matrix(|row, col| cofactor(matrix, col, row));

    // https://en.wikipedia.org/wiki/Minor_(linear_algebra)#Cofactor_expansion_of_the_determinant
    // Because we transposed the comatrix when we created it, we need to
    // calculate dot product of the first row of the matrix and first *column*
    // of the comatrix.
    let det = dot_product_with_column(&matrix[0], &comatrix_transposed, 0);
    if det.is_zero() {
        Err(crate::Error::DegenerateMatrix)
    } else {
        Ok((comatrix_transposed, det))
    }
}

#[test]
fn test_inverse_floats() {
    assert_eq!(
        Ok([
            [3.240813, -1.5373088, -0.49858665],
            [-0.96924317, 1.8759665, 0.041555062],
            [0.05563835, -0.20400737, 1.0571296]
        ]),
        inversed_copy(&crate::test::M_F32)
    );
}

#[test]
fn test_inverses_ratio() {
    assert_eq!(Ok(crate::test::M_INV_Q), inversed_copy(&crate::test::M_Q),);
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
    for<'x> &'x K: num_traits::NumOps<&'x K, K>,
{
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

/// Calculates a dot product of a 3-element vectors and a column in a matrix.
fn dot_product_with_column<K>(a: &[K; 3], b: &Matrix<K>, col: usize) -> K
where
    K: Scalar,
    for<'x> &'x K: num_traits::NumOps<&'x K, K>,
{
    &a[0] * &b[0][col] + &a[1] * &b[1][col] + &a[2] * &b[2][col]
}

#[test]
fn test_dot_product() {
    assert_eq!(4 + 10 + 18, dot_product(&[1i32, 2, 3], &[4i32, 5, 6]));
    let m = [[0, 4, 0], [0, 5, 0], [0, 6, 0]];
    assert_eq!(4 + 10 + 18, dot_product_with_column(&[1i32, 2, 3], &m, 1));
}

/// Returns a cofactor of a 3✕3 matrix M, i.e. C_{row,col}.
fn cofactor<K: Scalar>(matrix: &Matrix<K>, row: usize, col: usize) -> K
where
    for<'x> &'x K: num_traits::NumOps<&'x K, K>,
{
    let rr = ((row == 0) as usize, 2 - (row == 2) as usize);
    let cc = ((col == 0) as usize, 2 - (col == 2) as usize);
    let ad = &matrix[rr.0][cc.0] * &matrix[rr.1][cc.1];
    let bc = &matrix[rr.1][cc.0] * &matrix[rr.0][cc.1];
    if (row ^ col) & 1 == 0 {
        ad - bc
    } else {
        bc - ad
    }
}

#[test]
fn test_cofactor() {
    let matrix = [[1, 4, 7], [3, 0, 5], [-1, 9, 11]];
    let got = make_matrix(|row, col| cofactor(&matrix, row, col));
    assert_eq!([[-45, -38, 27], [19, 18, -13], [20, 16, -12]], got);
}


#[cfg(test)]
fn calc_diff(want: Matrix<crate::test::Rational>, got: Matrix<f32>) -> f64 {
    let mut err = num::BigRational::default();
    for (want, got) in want.into_iter().zip(got.into_iter()) {
        for (want, got) in want.into_iter().zip(got.into_iter()) {
            let want = want.into_raw();
            let want = num::BigRational::new_raw(want.0.into(), want.1.into());
            let got =  num::BigRational::from_float(got).unwrap();
            let e = want - got;
            err += &e * &e;
        }
    }
    num_traits::ToPrimitive::to_f64(&err).unwrap()
}

#[test]
#[ignore = "Expected to fail. The calculated error is the interesting part."]
fn test_accuracy() {
    let want = calculate(
        &crate::test::white(&crate::test::ratio).to_xyz(),
        &crate::test::primaries(&crate::test::ratio),
    )
    .unwrap();
    let got = calculate(
        &crate::test::white(&crate::test::float).to_xyz(),
        &crate::test::primaries(&crate::test::float),
    )
    .unwrap();
    let m_err = calc_diff(want, got);

    let want = inversed_copy(&want).unwrap();
    let got = inversed_copy(&got).unwrap();
    let inv_err = calc_diff(want, got);

    panic!("M   error: {:e}\nM⁻¹ error: {:e}", m_err, inv_err);
}
