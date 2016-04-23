library(testthat)
library(scran)

if (.Platform$OS.type!="windows") {
test_check("scran")
}
