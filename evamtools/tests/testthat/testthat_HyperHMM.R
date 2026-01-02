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

test_that("Check matrix dimensions",{
    set.seed(123) # Cualquier número sirve
    dataf <- data.frame(
            a = sample(c(0, 1), 10, replace = TRUE),
            b = sample(c(0, 1), 10, replace = TRUE),
            c = sample(c(0,1 ), 10, replace = TRUE))
    num_features <- ncol(dataf)
    expected_dim <- 2^num_features

    r_evam <- evam(dataf, method="HyperHMM")

    expect_equal(nrow(r_evam$HyperHMM_trans_mat),expected_dim)
    expect_equal(ncol(r_evam$HyperHMM_trans_mat),expected_dim)
})

test_that("HyperHMM transition matrix check", {
    set.seed(123)
    m <- data.frame(
    a = sample(c(0, 1), 10, replace = TRUE),
    b = sample(c(0, 1), 10, replace = TRUE),
    c = sample(c(0,1 ), 10, replace = TRUE)
    )
    
    r <- evam(m, method = "HyperHMM")
    
    # Accedemos a la matriz de transición (trans_mat o transition_rate_matrix)
    t_mat <- r$trans_matrix

    expect_s4_class(t_mat, "dgCMatrix") # Verifica que es una dgCMatrix

    states <- c("WT", "a", "b", "c", "a, b", "a, c", "b, c", "a, b, c")
    exp_mat <- matrix(0, nrow = 8, ncol = 8, dimnames = list(states, states))
    
    exp_mat["WT", "a"] <- 0.3364114
    exp_mat["WT", "b"] <- 0.4957018
    exp_mat["WT", "c"] <- 0.1678868

    exp_mat["a", "a, c"] <- 1.0000000
    exp_mat["a", "a, b"] <- 5.795146e-11 
    exp_mat["b", "b, c"] <- 1.0000000
    exp_mat["b", "a, b"] <- 9.166581e-10
    exp_mat["c", "a, c"] <- 0.8961278
    exp_mat["c", "b, c"] <- 0.1038722
    
    exp_mat["a, b", "a, b, c"] <- 1
    exp_mat["a, c", "a, b, c"] <- 1
    exp_mat["b, c", "a, b, c"] <- 1
    
    expect_equal(t_mat, exp_mat)
    
    expect_equal(dim(t_mat), c(8, 8))
})