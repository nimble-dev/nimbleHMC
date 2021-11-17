
test_that('tempFunction works', {
    expect_equal(tempFunction(0), 1)
    expect_equal(tempFunction(1), 2)
    expect_equal(tempFunction(2), 4)
})
