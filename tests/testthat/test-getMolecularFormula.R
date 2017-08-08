library(chem)

test_that('Chemical Formulae - rearrange', {
  out1 <- getMolecularFormula('C2H5OH', output = 'formula')
  expect_equal(out1, 'C2H6O')
})

test_that('Chemical Formulae - default Hill notation', {
  out2 <- getMolecularFormula('C3H4Be', output = 'formula')
  expect_equal(out2, 'C3H4Be')
})

test_that('Chemical Formulae - alphabetical', {
  out3 <- getMolecularFormula('C3H4Be', output = 'formula', order = 'alpha')
  expect_equal(out3, 'BeC3H4')
})

test_that('Chemical Formulae - empty', {
  out4A <- getMolecularFormula()
  out4B <- getMolecularFormula(NULL)
  out4C <- getMolecularFormula(NA)
  out4D <- getMolecularFormula('.')
  out4E <- getMolecularFormula('no formula')
  expect_equal(out4A, 'no formula')
  expect_equal(out4B, 'no formula')
  expect_equal(out4C, 'no formula')
  expect_equal(out4D, 'no formula')
  expect_equal(out4E, 'no formula')
})

test_that('Chemical Formulae - polymer', {
  out5A <- getMolecularFormula('(C3H4)n')
  out5B <- getMolecularFormula('polymer')
  expect_equal(out5A, 'polymer')
  expect_equal(out5B, 'polymer')
})

test_that('Chemical Formulae - R group', {
  out6A <- getMolecularFormula('C3H7R')
  out6B <- getMolecularFormula('R group')
  expect_equal(out6A, 'R group')
  expect_equal(out6B, 'R group')
})

test_that('Chemical Formulae - parentheses', {
  out7A <- getMolecularFormula('(CH3)3CH', output = 'formula')
  out7B <- getMolecularFormula('((CH3)3CH)2', output = 'formula')
  expect_equal(out7A, 'C4H10')
  expect_equal(out7B, 'C8H20')
})

test_that('Chemical Formulae - solvent', {
  out8A <- getMolecularFormula('(CH3)3CH.H2O', output = 'formula')
  out8B <- getMolecularFormula('(CH3)3CH.2H2O', output = 'formula')
  expect_equal(out8A, 'C4H12O')
  expect_equal(out8B, 'C4H14O2')
})

test_that('Chemical Formulae - table', {
  out9 <- getMolecularFormula('C2H5OH')
  expect_that(out9, is_a('data.frame'))
})
