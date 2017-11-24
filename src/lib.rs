extern crate num_traits;

use std::fmt::{self, Debug, Display, Formatter};
use std::cmp::Ordering;
use std::ops::{Add, Sub, Mul, Div};
use std::str::FromStr;

use self::num_traits::{Float, One, Zero};

/// A struct that represents positive floats by their natural logarithm.
///
/// # Examples
///
/// ```
/// # extern crate num_traits;
/// # extern crate log_domain;
/// # fn main() {
/// use log_domain::LogDomain;
/// use num_traits::{One, Zero};
///
/// match (LogDomain::new(0.5), LogDomain::new(0.25), LogDomain::new(0.75)) {
///     (Ok(x), Ok(y), Ok(z)) => {
///         assert_eq!(x.value(), 0.5);            // 0.5         = 0.5
///         assert_eq!(x.ln(), f64::ln(0.5));      // ln(0.5)     = ln(0.5)
///
///         // Operations `+`, `-`, `*`, `/`, and `pow`
///         assert_eq!(x + y, z);                  // 0.5  + 0.25 = 0.75
///         assert_eq!(z - x, y);                  // 0.75 - 0.5  = 0.25
///         assert_eq!(x * x, y);                  // 0.5  ⋅ 0.5  = 0.25
///         assert_eq!(y / x, x);                  // 0.25 / 0.5  = 0.5
///         assert_eq!(x.pow(2.0), y);             // 0.5²        = 0.25
///         assert_eq!(y.pow(1.0 / 2.0), x);       // √0.25       = 0.5
///
///         // Neutral elements `LogDomain::zero()` and `LogDomain::one()`
///         assert_eq!(z + LogDomain::zero(), z);  // 0.75 + 0    = 0.75
///         assert_eq!(z - z, LogDomain::zero());  // 0.75 - 0.75 = 0
///         assert_eq!(z * LogDomain::one(), z);   // 0.75 * 1    = 0.75
///         assert_eq!(z / z, LogDomain::one());   // 0.75 / 0.75 = 1
///         assert_eq!(z * LogDomain::zero(), LogDomain::zero());
///                                                // 0.75 * 0    = 0
///         assert_eq!(z.pow(0.0), LogDomain::one());
///                                                // 0.75⁰       = 1
///
///         // Comparison
///         assert!(z > y);                        // 0.75 > 0.25
///         assert!(y < z);                        // 0.25 < 0.75
///     },
///     _ => panic!(),
/// }
/// # }
/// ```
#[derive(PartialOrd, Debug, Clone, Copy)]
pub struct LogDomain<F: Float>(F);

impl<F: Float + Debug> LogDomain<F> {
    /// Creates a new `LogDomain` from a given value in the interval [0,∞).
    pub fn new(value: F) -> Result<Self, String> {
        if F::zero() <= value {
            Ok(LogDomain::new_unchecked(value))
        } else {
            Err(format!("{:?} is not in the interval [0,∞).", value))
        }
    }
}

impl<F: Float> LogDomain<F> {
    /// Logarithm of the value that is represented by the given `LogDomain`.
    pub fn ln(&self) -> F {
        match *self {
            LogDomain(value) => value,
        }
    }

    /// Same as `new`, but without bounds check.
    fn new_unchecked(value: F) -> Self {
        LogDomain(value.ln())
    }

    /// Value that is represented by the given `LogDomain`.
    pub fn value(&self) -> F {
        self.ln().exp()
    }

    /// Raise the represented probability to the power of a `Float` value.
    pub fn pow(&self, power: F) -> Self {
        match *self {
            LogDomain(value) => LogDomain(value * power),
        }
    }
}

/// An `impl` of `Ord` that defines
///  `Exp(NaN) = Exp(NaN)`,
///  `Exp(NaN) < Exp(y)`, and
///  `Exp(x) < Exp(y)` for `x < y`.
impl<F: Float> Ord for LogDomain<F> {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.partial_cmp(other) {
            Some(ordering) => ordering,
            None => {
                if self.ln().is_nan() {
                    if other.ln().is_nan() {
                        Ordering::Equal
                    } else {
                        Ordering::Greater
                    }
                } else {
                    Ordering::Less
                }
            }
        }
    }
}

/// An `impl` of `PartialEq` that defines
///  `Exp(NaN) = Exp(NaN)` and
///  `Exp(x) = Exp(y)` for `x - y < F::EPSILON`.
impl<F: Float> PartialEq for LogDomain<F> {
    fn eq(&self, other: &Self) -> bool {
        if self.ln().is_nan() {
            other.ln().is_nan()
        } else if other.ln().is_nan() {
            false
        } else {
            self.ln() == other.ln() || (self.ln() - other.ln()).abs() <= F::epsilon()
        }
    }
}

impl<F: Float> Eq for LogDomain<F> {}

/// An `impl` of `Add` that uses only two applications of transcendental functions
/// (`exp` and `ln_1p`) to increase precision.
impl<F: Float> Add for LogDomain<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let (a, b) = (self.ln(), other.ln());

        let (x, y) = if a > b {
            (a, b)
        } else {
            (b, a)
        };

        // Derivation of the formula:
        // Let x ≥ y.  Then
        // LogDomain(x) + LogDomain(y)
        // = exp(x) + exp(y)
        // = LogDomain( ln(exp(x) + exp(y)) )
        // = LogDomain( ln(exp(x) + exp(x) ⋅ exp(y) / exp(x)) )
        // = LogDomain( ln(exp(x) + exp(x) ⋅ exp(y - x)) )
        // = LogDomain( ln(exp(x) ⋅ (1 + exp(y - x))) )
        // = LogDomain( x + ln(1 + exp(y - x)) )
        // = LogDomain( x + ln_1p(exp(y - x)) )
        LogDomain(x + (y - x).exp().ln_1p())
    }
}

/// An `impl` of `Sub` that uses only two applications of transcendental functions
/// (`exp_m1` and `ln`) to increase precision.
impl<F: Float + Debug> Sub for LogDomain<F> {
    type Output = Self ;

    fn sub(self, other: Self) -> Self {
        match (self.ln(), other.ln()) {
            // Derivation of the formula:
            // Let x > y. Then
            // LogDomain(x) - LogDomain(y)
            // = exp(x) - exp(y)
            // = LogDomain( ln(exp(x) - exp(y)) )
            // = LogDomain( ln(exp(x) - exp(x) ⋅ exp(y) / exp(x)) )
            // = LogDomain( ln(exp(x) - exp(x) ⋅ exp(y - x)) )
            // = LogDomain( x + ln(1 - exp(y - x)) )
            // = LogDomain( x + ln(- exp_m1(y - x)) )
            (x, y) if x > y  => LogDomain(x + (-(y - x).exp_m1()).ln()),
            (x, y) if x == y => LogDomain::zero(),
            (x, y) if x <  y => panic!("exp({:?}) - exp({:?}) is less than zero", x, y),
            _                => unreachable!(),
        }
    }
}

impl<F: Float> Mul for LogDomain<F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        LogDomain(self.ln().add(other.ln()))
    }
}

impl<F: Float> Div for LogDomain<F> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        if !other.is_zero() {
            LogDomain(self.ln().sub(other.ln()))
        } else {
            panic!("division by zero")
        }
    }
}

impl<F: Float> Zero for LogDomain<F> {
    fn zero() -> Self {
        LogDomain(F::neg_infinity())
    }

    fn is_zero(&self) -> bool {
        self.ln() == F::neg_infinity()
    }
}

impl<F: Float> One for LogDomain<F> {
    fn one() -> Self {
        LogDomain(F::zero())
    }
}

impl<F: Debug + Float + FromStr<Err=E>, E: ToString> FromStr for LogDomain<F> {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.parse() {
            Ok(p) => LogDomain::new(p),
            Err(e) => Err(e.to_string()),
        }
    }
}

impl<F: Float + Display> Display for LogDomain<F> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{}", self.value())
    }
}
