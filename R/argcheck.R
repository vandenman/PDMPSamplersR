validate_type <- function(
    x,
    type = c("double", "integer", "character", "logical"),
    n = NULL,
    dims = NULL,
    allow_missing = FALSE,
    allow_inf = FALSE,
    positive = FALSE,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()
) {
  rlang::check_required(x, arg = arg, call = call)

  # Match argument type
  type <- match.arg(type, several.ok = TRUE)

  # Actual type
  actual_type <- typeof(x)
  if (!(actual_type %in% type)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be of type {.val {type}}.",
        x = "Got {.val {actual_type}}."
      ),
      call = call
    )
  }

  # Missingness
  if (!allow_missing && anyNA(x)) {
    cli::cli_abort("{.arg {arg}} contains missing values but they are disallowed.", call = call)
  }

  # Infinite values (only makes sense for numeric types)
  if (!allow_inf && actual_type %in% c("double", "integer") && any(is.infinite(x))) {
    cli::cli_abort("{.arg {arg}} contains infinite values but they are disallowed.", call = call)
  }

  # Positive values (only for numeric types)
  if (positive && actual_type %in% c("double", "integer")) {
    if (any(x <= 0, na.rm = TRUE)) {
      cli::cli_abort("{.arg {arg}} must contain only positive values.", call = call)
    }
  } else if (positive && !(actual_type %in% c("double", "integer"))) {
    cli::cli_abort("{.arg {arg}} positivity check is only valid for numeric types.", call = call)
  }

  # Length check (only for vectors)
  if (!is.null(n) && is.null(dim(x))) {
    if (length(x) != n) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must be length {n}.",
          x = "Got length {length(x)}."
        ),
        call = call
      )
    }
  }

  # Dimension check (for arrays/matrices)
  if (!is.null(dims)) {
    if (is.null(dim(x))) {
      cli::cli_abort(
        "{.arg {arg}} must be an array/matrix with dimensions {dims}, but has no dimensions.",
        call = call
      )
    }
    actual_dims <- dim(x)

    if (length(dims) != length(actual_dims)) {
      cli::cli_abort(
        c(
          "Invalid dimensions for {.arg {arg}}.",
          i = "Expected {.val {dims}} (length {length(dims)}).",
          x = "Got {.val {actual_dims}} (length {length(actual_dims)})."
        ),
        call = call
      )
    }

    mismatch <- !is.na(dims) & (dims != actual_dims)
    if (any(mismatch)) {
      cli::cli_abort(
        c(
          "Invalid dimensions for {.arg {arg}}.",
          i = "Expected {.val {dims}} (NA means 'any').",
          x = "Got {.val {actual_dims}}."
        ),
        call = call
      )
    }
  }

  invisible(TRUE)
}


cast_integer <- function(x, n = NULL) {
  if (rlang::is_integerish(x, n = n))
    x <- as.integer(x)
  x
}