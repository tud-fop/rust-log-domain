# rust-log-domain
Logarithmic representation of floats for rust.

## License

BSD-3-clause

## Documentation
* see `src/lib.rs`
* generate HTML-documentation with `cargo doc`

## Use
* Include the crate in your rust project by adding
  ```
  log_domain = { git = "https://github.com/tud-fop/rust-log-domain.git" }
  ```
  to the `[dependencies]` in your `Cargo.toml`.
* The crate contains a `struct LogDomain<F>` and appropriate implementations for common arithmetic operations.  The struct is a newtype to allow zero-cost abstraction.
