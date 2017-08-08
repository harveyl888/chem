library(chem)

test_that('Exact Mass', {
  out1A <- exactMass('H2O')  # water
  out1B <- exactMass('C27H46O')  # cholersterol
  expect_equal(out1A, 18.0106, tolerance = 1E-4)
  expect_equal(out1B, 386.35486, tolerance = 1E-4)
})

test_that('Exact Mass - missing atoms', {
  out2 <- exactMass('C2H5Q3')
  expect_equal(out2, 'atoms not recognized: Q')
})
