extern crate num_traits;

use std::{fmt, f64};
use std::cmp::Ordering;
use std::ops::{Add, Sub, Mul, Div};
use std::str::FromStr;

use self::num_traits::{One, Zero};

/// A struct that represents probabilities.
/// A probability is internally represented as its natural logarithm.
/// Probabilities greater than 1 are allowed during calculations.
///
/// # Examples
///
/// ```
/// # extern crate num_traits;
/// # extern crate log_prob;
/// # fn main() {
/// use log_prob::LogProb;
/// use num_traits::{One, Zero};
///
/// match (LogProb::new(0.5), LogProb::new(0.25), LogProb::new(0.75)) {
///     (Ok(x), Ok(y), Ok(z)) => {
///         assert_eq!(x.probability(), 0.5);    // 0.5         = 0.5
///         assert_eq!(x.ln(), f64::ln(0.5));    // ln(0.5)     = ln(0.5)
///
///         // Operations `+`, `-`, `*`, and `/`
///         assert_eq!(x + y, z);                // 0.5  + 0.25 = 0.75
///         assert_eq!(z - x, y);                // 0.75 - 0.5  = 0.25
///         assert_eq!(x * x, y);                // 0.5  ⋅ 0.5  = 0.25
///         assert_eq!(y / x, x);                // 0.25 / 0.5  = 0.5
///
///         // Neutral elements `LogProb::zero()` and `LogProb::one()`
///         assert_eq!(z + LogProb::zero(), z);  // 0.75 + 0    = 0.75
///         assert_eq!(z - z, LogProb::zero());  // 0.75 - 0.75 = 0
///         assert_eq!(z * LogProb::one(), z);   // 0.75 * 1    = 0.75
///         assert_eq!(z / z, LogProb::one());   // 0.75 / 0.75 = 1
///         assert_eq!(z * LogProb::zero(), LogProb::zero());
///                                              // 0.75 * 0    = 0
///
///         // Comparison
///         assert!(z > y);                      // 0.75 > 0.25
///         assert!(y < z);                      // 0.25 < 0.75
///     },
///     _ => panic!(),
/// }
/// # }
/// ```
#[derive(PartialOrd, Debug, Clone, Copy)]
pub struct LogProb(f64);

impl LogProb {
    /// Creates a new `LogProb` from a given value in the interval [0,∞).
    pub fn new(value: f64) -> Result<Self, String> {
        if 0.0 <= value {
            Ok(LogProb::new_unchecked(value))
        } else {
            Err(format!("{} is not a probability, i.e. not in the interval [0,∞).", value))
        }
    }

    /// Logarithm of the probability that is represented by the given `LogProb`.
    pub fn ln(&self) -> f64 {
        match self {
            &LogProb(value) => value,
        }
    }

    /// Same as `new`, but without bounds check.
    pub fn new_unchecked(value: f64) -> Self {
        LogProb(value.ln())
    }

    /// Probability that is represented by the given `LogProb`.
    pub fn probability(&self) -> f64 {
        self.ln().exp()
    }
}

/// An `impl` of `Ord` that defines `Exp(NaN) = Exp(NaN)`, `Exp(NaN) < Exp(y)`, and `Exp(x) < Exp(y)` for `x < y`.
impl Ord for LogProb {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.partial_cmp(&other) {
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

/// An `impl` of `PartialEq` that defines `Exp(NaN) = Exp(NaN)` and `Exp(x) = Exp(y)` for `x - y < f64::EPSILON`.
impl PartialEq for LogProb {
    fn eq(&self, other: &Self) -> bool {
        if self.ln().is_nan() {
            if other.ln().is_nan() { true } else { false }
        } else if other.ln().is_nan() {
            false
        } else {
            self.ln() == other.ln() || (self.ln() - other.ln()).abs() <= f64::EPSILON
        }
    }
}

impl Eq for LogProb {}

/// An `impl` of `Add` that uses only two applications of transcendental functions (`exp` and `ln_1p`) to increase precision.
impl Add for LogProb {
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
        // LogProb(x) + LogProb(y)
        // = exp(x) + exp(y)
        // = LogProb( ln(exp(x) + exp(y)) )
        // = LogProb( ln(exp(x) + exp(x) ⋅ exp(y) / exp(x)) )
        // = LogProb( ln(exp(x) + exp(x) ⋅ exp(y - x)) )
        // = LogProb( ln(exp(x) ⋅ (1 + exp(y - x))) )
        // = LogProb( x + ln(1 + exp(y - x)) )
        // = LogProb( x + ln_1p(exp(y - x)) )
        LogProb(x + (y - x).exp().ln_1p())
    }
}

/// An `impl` of `Sub` that uses only two applications of transcendental functions (`exp_m1` and `ln`) to increase precision.
impl Sub for LogProb {
    type Output = Self ;

    fn sub(self, other: Self) -> Self {
        match (self.ln(), other.ln()) {
            // Derivation of the formula:
            // Let x > y. Then
            // LogProb(x) - LogProb(y)
            // = exp(x) - exp(y)
            // = LogProb( ln(exp(x) - exp(y)) )
            // = LogProb( ln(exp(x) - exp(x) ⋅ exp(y) / exp(x)) )
            // = LogProb( ln(exp(x) - exp(x) ⋅ exp(y - x)) )
            // = LogProb( x + ln(1 - exp(y - x)) )
            // = LogProb( x + ln(- exp_m1(y - x)) )
            (x, y) if x > y  => LogProb(x + (-(y - x).exp_m1()).ln()),
            (x, y) if x == y => LogProb::zero(),
            (x, y) if x <  y => panic!("exp({}) - exp({}) is less than zero", x, y),
            _                => unreachable!(),
        }
    }
}

impl Mul for LogProb {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        LogProb(self.ln().add(other.ln()))
    }
}

impl Div for LogProb {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        if !other.is_zero() {
            LogProb(self.ln().sub(other.ln()))
        } else {
            panic!("division by zero")
        }
    }
}

impl Zero for LogProb {
    fn zero() -> LogProb {
        LogProb(f64::NEG_INFINITY)
    }

    fn is_zero(&self) -> bool {
        self.ln() == f64::NEG_INFINITY
    }
}

impl One for LogProb {
    fn one() -> LogProb {
        LogProb(0.0)
    }
}

impl<T: Into<f64>> From<T> for LogProb {
    fn from(value: T) -> Self {
        let the_f64 = value.into();
        if 0.0 <= the_f64 {
            LogProb::new_unchecked(the_f64)
        } else {
            panic!("{} is not a probability, i.e. not in the interval [0,∞).", the_f64)
        }
    }
}

impl FromStr for LogProb {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.parse() {
            Ok(p) => LogProb::new(p),
            Err(e) => Err(e.to_string())
        }
    }
}

impl fmt::Display for LogProb {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.probability())
    }
}
