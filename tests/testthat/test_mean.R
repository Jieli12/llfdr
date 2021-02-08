context("Test Mean Estimator")

test_that("Mean Function", {
    data(Y)
    data(u)
    data(h1_lc)
    ylc <- YResidualH1LC(Y = Y, u = u, h = h1_lc$minimum)
    expect_identical(class(ylc)[1], "matrix")
})
