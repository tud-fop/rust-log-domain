# rust-log-domain
Logarithmic representation of floats for rust.

## License

BSD-3-clause

## Documentation
* https://docs.rs/log_domain/0.4.0/log_domain/
* or generate HTML-documentation with `cargo doc`

## Use
* Include the crate in your rust project by adding
  ```
  log_domain = "0.4.0"
  ```
  to the `[dependencies]` in your `Cargo.toml`.
* The crate contains a `struct LogDomain<F>` and appropriate implementations for common arithmetic operations.  The struct is a newtype to allow zero-cost abstraction.
