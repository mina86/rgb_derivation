# RGB colour system derivation routines

Functions for deriving RGB→XYZ and XYZ→RGB conversion matrices for
given RGB colour system (such as sRGB colour space).  The calculations
are performed from the definition of such system provided in the form
of chromacicities of the reference white point and the red, green and
blue primaries.  Alternatively, constructions from XYZ coordinates of
primaries is also available.

The crate supports arithmetic with rational and big integer types such
that the calculations can be performed without any loss of precision
if desired.  So long as a type implements the four basic arithmetic
operations, it can be used with this library.  For example, `f32`,
`num::Rational64` and `num::BigRational` can all be used.


## Usage

Using this package with Cargo projects is as simple as adding a single
dependency:

```toml
[dependencies]
rgb_derivation = "0.2"
```

With that dependency in place, it’s now simple to write an application
which converts an sRGB colour into other colour spaces:


```rust
type Scalar = num::rational::Ratio<i128>;
type Chromaticity = rgb_derivation::Chromaticity<Scalar>;

fn chromaticity(x: (i128, i128), y: (i128, i128)) -> Chromaticity {
    Chromaticity::new(Scalar::new(x.0, x.1), Scalar::new(y.0, y.1)).unwrap()
}

fn print_vector(header: &str, vector: &[Scalar; 3]) {
    print!("{}: [", header);
    for (idx, value) in vector.iter().enumerate() {
        print!("{} {} / {}",
               if idx == 0 { "" } else { "," },
               value.numer(), value.denom());
    }
    println!(" ]");
}

fn print_matrix(header: &str, matrix: &[[Scalar; 3]; 3]) {
    static OPEN: [char; 3] = ['⎡', '⎢', '⎣'];
    static CLOSE: [char; 3] = ['⎤', '⎥', '⎦'];

    fn make_array<T>(f: impl Fn(usize) -> T) -> [T; 3] { [f(0), f(1), f(2)] }

    let formatted = make_array(|row| make_array(|col| {
        let value = &matrix[row][col];
        (format!("{}", value.numer()), format!("{}", value.denom()))
    }));
    let lengths = make_array(|col| (
        formatted.iter().map(|row| row[col].0.len()).max().unwrap(),
        formatted.iter().map(|row| row[col].1.len()).max().unwrap(),
    ));

    let indent = header.chars().count();
    for (idx, row) in formatted.iter().enumerate() {
        if idx == 1 {
            print!("{}: {}", header, OPEN[idx]);
        } else {
            print!("{:indent$}  {}", "", OPEN[idx], indent = indent);
        }
        for (idx, value) in row.iter().enumerate() {
            print!("{comma} {numer:>numer_len$} / {denom:>denom_len$}",
                   comma = if idx == 0 { "" } else { "," },
                   numer = value.0, numer_len = lengths[idx].0,
                   denom = value.1, denom_len = lengths[idx].1);
        }
        println!(" {}", CLOSE[idx]);
    }
}

fn main() {
    let white_xy = chromaticity((312713, 1000000), (329016, 1000000));
    let primaries_xy = [
        chromaticity((64, 100), (33, 100)),
        chromaticity((30, 100), (60, 100)),
        chromaticity((15, 100), (6, 100)),
    ];

    let white_xyz = white_xy.to_xyz();
    let matrix = rgb_derivation::matrix::calculate(
        &white_xyz, &primaries_xy).unwrap();
    let inverse = rgb_derivation::matrix::inversed_copy(&matrix).unwrap();
    let primaries_xyz = rgb_derivation::matrix::transposed_copy(&matrix);

    print_vector("sRGB white point (D65)", &white_xyz);
    print_matrix("sRGB primaries", &primaries_xyz);
    print_matrix("sRGB→XYZ", &matrix);
    print_matrix("XYZ→sRGB", &inverse);
}
```
