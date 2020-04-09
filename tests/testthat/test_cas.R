library(chem)

test_that('CAS Check', {
  out1A <- cas_check('64-17-5')  # ethanol
  out1B <- cas_check('64-17-6')  # incorrect cas
  expect_true(out1A)
  expect_false(out1B)
})

test_that('Retrieve smiles from CAS', {
  out1A <- cas_to_smiles('64-17-5')  # ethanol
  expect_equal(out1A, 'CCO')
  expect_warning(cas_to_smiles('64-17-6'))  # incorrect cas
})
