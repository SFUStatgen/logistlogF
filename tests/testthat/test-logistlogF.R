test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
test_that("logistlogF works on DES data",{
  DES$fmatched <- factor(DES$matched.set)
  fit <- logistlogF(case~fmatched+DES+matern.smoke,dat=DES,m=2)
  ans <- c(-3.22357637,-0.12714946,-0.07236834,-0.12714946,
           0.82809575,-0.12714946,-0.04467054,-0.12714946,
           4.70333417,0.50460302)
  names(ans) <- c("Intercept","fmatched2","fmatched3","fmatched4","fmatched5",
                  "fmatched6","fmatched7","fmatched8","DES","matern.smoke")
  expect_equal(fit$coefficients,ans)
})
