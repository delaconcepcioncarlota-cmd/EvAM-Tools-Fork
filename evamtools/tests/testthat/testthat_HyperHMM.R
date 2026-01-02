set.seed(123) # Cualquier número sirve
datos <- data.frame(
  a = sample(c(0, 1), 10, replace = TRUE),
  b = sample(c(0, 1), 10, replace = TRUE),
  c = sample(c(0,1 ), 10, replace = TRUE)
)

#Check of translation from binary to int
test_that("Transforming for binary to integer",{
    set.seed(123) # Cualquier número sirve
    dataf <- data.frame(
            a = sample(c(0, 1), 10, replace = TRUE),
            b = sample(c(0, 1), 10, replace = TRUE),
            c = sample(c(0,1 ), 10, replace = TRUE))
    datam <- as.matrix(dataf)
    feature_labels <- colnames(dataf)
    num_features <- ncol(dataf)

    r <- evam(datam, methods="HyperHMM")

    all_ids <- 0:(2^(num_features)-1)
    decoded_states <- vapply(all_ids, translate_state,
                            character(1),
                            genes = feature_labels)
    expect_equal(translate_state(0,genes=feature_labels),"WT")
    expect_equal(translate_state(1,genes=feature_labels),"c")
    expect_equal(translate_state(2,genes=feature_labels),"b")
    expect_equal(translate_state(3,genes=feature_labels),"b, c")
    expect_equal(translate_state(4,genes=feature_labels),"a")
    expect_equal(translate_state(5,genes=feature_labels),"a, c")
    expect_equal(translate_state(6,genes=feature_labels),"a, b")
    expect_equal(translate_state(7,genes=feature_labels),"a, b, c")
    expect_equal(translate_state(8,genes=feature_labels),"")
    }
)

