# RGB colour system derivation routines

Functions for deriving RGB→XYZ and XYZ→RGB conversion matrices for given RGB
colour system (such as sRGB colour space).  The calculations are performed from
the definition of such system provided in the form of chromacicities of the
reference white point and the red, green and blue primaries.  Alternatively,
constructions from XYZ coordinates of primaries is also available.

The crate supports calculations using any numeric type which handles the four
basic arithmetic operations.  In particular, rational types such as
`num::rational::Ratio<i128>` or `num::BigRational` can be used to avoid any loss
of precision while performing the calculations.  (Note that `num::Rational64`
and especially `num::Rational32` may lead to overflows thus those types are not
recommended unless overflow checking is enabled).

Of course the calculations can also be performed with `f32` or `f64` primitive
types.


## Usage

Using this package with Cargo projects is as simple as adding a single
dependency:

```toml
[dependencies]
rgb_derivation = "0.2.2"
```

With that dependency in place, it’s now simple to write an application
which converts an sRGB colour into other colour spaces:

```rust
type Scalar = num::rational::Ratio<i128>;
type Chromaticity = rgb_derivation::Chromaticity<Scalar>;

fn chromaticity(x: (i128, i128), y: (i128, i128)) -> Chromaticity {
    let x = Scalar::new(x.0, x.1);
    let y = Scalar::new(y.0, y.1);
    Chromaticity::new(x, y).unwrap()
}

fn print_vector(header: &str, vector: &[Scalar; 3]) {
    print!("{}: [", header);
    for (idx, value) in vector.iter().enumerate() {
        let sep = if idx == 0 { "" } else { "," };
        print!("{sep} {} / {}", value.numer(), value.denom());
    }
    println!(" ]");
}

fn print_matrix(header: &str, matrix: &[[Scalar; 3]; 3]) {
    static OPEN: [char; 3] = ['⎡', '⎢', '⎣'];
    static CLOSE: [char; 3] = ['⎤', '⎥', '⎦'];

    let formatted = matrix.map(|row| {
        row.map(|v| (v.numer().to_string(), v.denom().to_string()))
    });
    let lengths = [0, 1, 2].map(|col| (
        formatted.iter().map(|row| row[col].0.len()).max().unwrap(),
        formatted.iter().map(|row| row[col].1.len()).max().unwrap(),
    ));

    let padding = header.chars().count();
    for (idx, row) in formatted.iter().enumerate() {
        if idx == 1 {
            print!("{header}: {}", OPEN[idx]);
        } else {
            print!("{:padding$}  {}", "", OPEN[idx]);
        }
        let mut comma = "";
        for (value, lens) in row.iter().zip(lengths.iter()) {
            let (numer, denom) = value;
            let (numer_len, denom_len) = lens;
            print!("{comma} {numer:>numer_len$} / {denom:>denom_len$}");
            comma = ",";
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

    let white_xyz = white_xy.into_xyz();
    let (matrix, inverse) =
        rgb_derivation::calculate_pair(&white_xyz, &primaries_xy)
            .unwrap();
    let primaries_xyz =
        rgb_derivation::matrix::transposed_copy(&matrix);

    print_vector("sRGB white point (D65)", &white_xyz);
    print_matrix("sRGB primaries", &primaries_xyz);
    print_matrix("sRGB→XYZ", &matrix);
    print_matrix("XYZ→sRGB", &inverse);
}
```

Note: if you need matrices for the sRGB colour space, the [`srgb`
crate](https://crates.io/crates/srgb) provides them along with gamma functions
needed to properly handle sRGB.
