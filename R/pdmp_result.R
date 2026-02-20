# Constructor (internal)
new_pdmp_result <- function(chains, stats, d, n_chains) {
  structure(
    list(
      chains   = chains,
      stats    = stats,
      d        = d,
      n_chains = n_chains
    ),
    class = "pdmp_result"
  )
}

#' @export
print.pdmp_result <- function(x, ...) {
  cat(sprintf(
    "PDMP result: %d dimension%s, %d chain%s\n",
    x$d, if (x$d == 1L) "" else "s",
    x$n_chains, if (x$n_chains == 1L) "" else "s"
  ))
  invisible(x)
}
