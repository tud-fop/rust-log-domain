# rust-log-prob
Logarithmic probabilities for rust.

## License

BSD-3-clause

## Documentation
* see `src/lib.rs`
* generate HTML-documentation with `cargo doc`

## Use
* Include the crate in your rust project by adding
  ```
  log_prob = { git = "https://github.com/tud-fop/rust-log-prob.git" }
  ```
  to the `[dependencies]` in your `Cargo.toml`.
* The crate contains a `struct LogProb` and appropriate implementations for common arithmetic operations.  The struct is a newtype to allow zero-cost abstraction.
