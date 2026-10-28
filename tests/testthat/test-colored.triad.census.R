
#' check that non-zero elements of named vector match expected values
#'
#' @param actual result of a triad census
#' @param expected named vector of expected values
#' @examples
#' # will pass, without having to specify foo=0
#' expect_equal_census(c(foo=0, bar=2), c(bar=2))
#' \dontrun{
#' # will fail
#' expect_equal_census(c(foo=0, bar=2), c(bar=1))
#' expect_equal_census(c(foo=0, bar=2), c(bar=2, bam=0))
#' }
expect_equal_census <- function(actual, expected) {

  # homogenize names
  bad  <- setdiff(names(expected), names(actual))
  actual[bad] <- NA

  bad  <- setdiff(names(actual), names(expected))
  expected[bad] <- 0

  expected <- expected[names(actual)]

  # report problematic
  bad <- actual != expected | is.na(actual) | is.na(expected)
  if (any(bad)) {
    msg <- mapply(
      function(n, x, y) paste0(n, ': actual=', x, ', expected=', y),
      names(bad)[bad], actual[bad], expected[bad]
    )
    fail(sprintf("Triads with differing values:\n\t%s", paste(msg, collapse = "\n\t")))
  }

  succeed()
  invisible(actual)
}

test_that("T003-111", {
  triads <- colored.triad.census(matrix(0,3,3), c(1,1,1))
  expect_equal_census(triads, c("T003-111" = 1))

  triads <- colored.triad.census(matrix(0,3,3), c(1,1,1), directed = TRUE)
  expect_equal_census(triads, c("T003-111" = 1))

  triads <- colored.triad.census(matrix(0,4,4), c(2,2,2,2))
  expect_equal_census(triads, c("T003-111" = 4))
})

test_that("T003-123", {
  triads <- colored.triad.census(matrix(0,3,3), c(1,2,3))
  expect_equal_census(triads, c("T003-123" = 1))

  triads <- colored.triad.census(matrix(0,3,3), c(3,2,1))
  expect_equal_census(triads, c("T003-123" = 1))
})
