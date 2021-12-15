
test_that("Phenotype name generation works", {
  expect_equal(make_phenotype_name(c(0,1,2,-1),c("Marker1","Marker2","Marker3","Marker4")), "Marker1-Marker2+Marker3++")
})

test_that("Marker count works", {
  expect_equal(count_markers(c(0,1,2,-1)),3)
})
