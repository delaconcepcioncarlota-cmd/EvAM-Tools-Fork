t1 <- Sys.time()

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
    expect_equal(evamtools:::translate_state(0,genes=feature_labels),"WT")
    expect_equal(evamtools:::translate_state(1,genes=feature_labels),"c")
    expect_equal(evamtools:::translate_state(2,genes=feature_labels),"b")
    expect_equal(evamtools:::translate_state(3,genes=feature_labels),"b, c")
    expect_equal(evamtools::translate_state(4,genes=feature_labels),"a")
    expect_equal(evamtools:::translate_state(5,genes=feature_labels),"a, c")
    expect_equal(evamtools:::translate_state(6,genes=feature_labels),"a, b")
    expect_equal(evamtools:::translate_state(7,genes=feature_labels),"a, b, c")
    expect_equal(evamtools:::translate_state(8,genes=feature_labels),"")
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

    r_evam <- evam(dataf, methods="HyperHMM")

    expect_equal(nrow(r_evam$HyperHMM_trans_mat),expected_dim)
    expect_equal(ncol(r_evam$HyperHMM_trans_mat),expected_dim)
})

test_that("Check number of features",{
    set.seed(123)
    dataf <- data.frame(
            a = sample(c(0, 1), 10, replace = TRUE),
            b = sample(c(0, 1), 10, replace = TRUE),
            c = sample(c(0,1 ), 10, replace = TRUE))
    datam <- as.matrix(dataf)
    feature_labels <- colnames(dataf)
    num_features <- ncol(dataf)
    expected_dim <- 2^num_features

    r_evam <- evam(datam, methods="HyperHMM")
    r_hhmm <- HyperHMM::hyperhmm(datam)

    expect_equal(r_evam$HyperHMM_n_features$N, num_features)
    expect_equal(r_hhmm$L, num_features)
    expect_identical(r_evam$HyperHMM_n_features$N, r_hhmm$L)
})

test_that("HyperHMM transition matrix check", {
    set.seed(123)
    m <- data.frame(
    a = sample(c(0, 1), 10, replace = TRUE),
    b = sample(c(0, 1), 10, replace = TRUE),
    c = sample(c(0,1 ), 10, replace = TRUE))
    
    r <- evam(m, methods="HyperHMM")
    
    t_mat <- r$HyperHMM_trans_mat

    expect_true(inherits(t_mat, "dgCMatrix"))

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

    expect_equal(as.matrix(t_mat), exp_mat, tolerance = 1e-6)

    expect_equal(dim(t_mat), c(8, 8))
})

test_that("HyperHMM transition rate matrix check", {
    set.seed(123)

    m <- data.frame(
    a = sample(c(0, 1), 10, replace = TRUE),
    b = sample(c(0, 1), 10, replace = TRUE),
    c = sample(c(0,1 ), 10, replace = TRUE)
    )
    
    r <- evam(m, methods="HyperHMM")
    
    t_rate_mat <- r$HyperHMM_trans_rate_mat

    expect_true(inherits(t_rate_mat, "dgCMatrix"))

    states <- c("WT", "a", "b", "c", "a, b", "a, c", "b, c", "a, b, c")
    exp_rate_mat <- matrix(0, nrow = 8, ncol = 8, dimnames = list(states, states))

    exp_rate_mat["WT", "a"] <- 0.3364114
    exp_rate_mat["WT", "b"] <- 0.4957018
    exp_rate_mat["WT", "c"] <- 0.1678868
    
    exp_rate_mat["a", "a, c"] <- 0.3364114
    exp_rate_mat["a", "a, b"] <- 1.949553e-11
    exp_rate_mat["b", "b, c"] <- 0.49570178
    exp_rate_mat["b", "a, b"] <- 4.543890e-10
    exp_rate_mat["c", "a, c"] <- 0.1504480
    exp_rate_mat["c", "b, c"] <- 0.01743878
    
    exp_rate_mat["a, b", "a, b, c"] <- 4.738846e-10
    exp_rate_mat["a, c", "a, b, c"] <- 4.868594e-01
    exp_rate_mat["b, c", "a, b, c"] <- 5.131406e-01

    expect_equal(as.matrix(t_rate_mat), exp_rate_mat, tolerance = 1e-6)

    expect_equal(dim(t_rate_mat), c(8, 8))
})

test_that("Check initialstates argument on hyper_hmm_opts function",{
    m = matrix(c(1,0,0,0,0,
                   1,1,0,0,0,
                   1,1,1,0,0,
                   1,1,1,1,0,
                   1,1,1,1,1,
                   0,0,0,0,1,
                   0,0,0,1,1,
                   0,0,1,1,1,
                   0,1,1,1,1,
                   1,1,1,1,1), byrow=TRUE, ncol = 5)
    colnames(m)<-c("a","b","c","d","e")
    #Matriz de estados iniciales
    m_initial_states = matrix(c(0,0,0,0,0,
                   1,0,0,0,0,
                   1,1,0,0,0,
                   1,1,1,0,0,
                   1,1,1,1,0,
                   0,0,0,0,0,
                   0,0,0,0,1,
                   0,0,0,1,1,
                   0,0,1,1,1,
                   0,1,1,1,1), byrow=TRUE, ncol = 5)
    colnames(m_initial_states)<-c("a","b","c","d","e")

    feature_labels <- colnames(m)
    num_features <- ncol(m)

    #Opciones pasadas a evam como argumento (lista)
    opts_HyperHMM = list(initialstates=m_initial_states)

    r2<-evam(m, methods="HyperHMM", hyper_hmm_opts=opts_HyperHMM)
    t_mat2 <- r2$HyperHMM_trans_mat
    #Test para comprobar que genera una matriz de transición correctamente
    expect_s4_class(t_mat2, "dgCMatrix")
})

test_that("HyperHMM plotting basic execution", {
    set.seed(123)

    dataf <- data.frame(
    a = sample(c(0, 1), 10, replace = TRUE),
    b = sample(c(0, 1), 10, replace = TRUE),
    c = sample(c(0,1 ), 10, replace = TRUE)
    )
    r <- evam(as.matrix(dataf), methods = "HyperHMM")
    
    # Verificación de objetos gráficos
    expect_s3_class(plot_HyperHMM_bubbles(r$HyperHMM_primary_output), "ggplot")
    expect_s3_class(plot_HyperHMM_hypercube(r$HyperHMM_primary_output), "ggplot")
    expect_silent(suppressMessages(plot_HyperHMM_pfg(r$HyperHMM_primary_output, pfg.layout = "matrix")))
})

test_that("Check only correct input dimensions", {

    m <- matrix(c(1, 0, 0, 1, 1, 1), ncol = 2)

    colnames(m) <- c("a", "b")

    expect_error(evam(m, method = "HyperHMM"))
    expect_error(evam(c(1, 0, 1), method = "HyperHMM"))
})

test_that("HyperHMM plotting functions handle different thresholds", {
    set.seed(123)
    dataf <- data.frame(
        a = sample(c(0, 1), 10, replace = TRUE),
        b = sample(c(0, 1), 10, replace = TRUE),
        c = sample(c(0, 1), 10, replace = TRUE)
    )
    
    r <- evam(dataf, methods = "HyperHMM")

    expect_no_error(plot_HyperHMM_hypercube(r$HyperHMM_primary_output, threshold = 0.99))

    p_flux <- plot_HyperHMM_hypercube_flux(r$HyperHMM_primary_output, thresh = 0.01)
    expect_s3_class(p_flux, "ggplot")
})

test_that("Check evam function dealing with spelling errors in HyperHMM", {
    set.seed(123)
    dataf <- data.frame(
            a = sample(c(0, 1), 10, replace = TRUE),
            b = sample(c(0, 1), 10, replace = TRUE),
            c = sample(c(0,1 ), 10, replace = TRUE))
    datam <- as.matrix(dataf)

    try_correct=evam(datam, methods="HyperHMM")
    try_Hyperhmm=evam(datam, methods="Hyperhmm")
    try_HyperHmm=evam(datam, methods="hyperHmm")
    try_hyperhMm=evam(datam, methods="hyperhMm")
    try_hyperhmM=evam(datam, methods="hyperhmM")
    try_HyperHmm=evam(datam, methods="HyperHmm")
    try_HyperHMm=evam(datam, methods="HyperHMm")
    try_HyperhMM=evam(datam, methods="HyperhMM")
    try_HyperHmM=evam(datam, methods="HyperHmM")
    try_hyperHMM=evam(datam, methods="hyperHMM")

    expect_equal(try_correct$HyperHMM_trans_mat, try_Hyperhmm$HyperHMM_trans_mat)
    expect_equal(try_correct$HyperHMM_trans_mat, try_HyperHmm$HyperHMM_trans_mat)
    expect_equal(try_correct$HyperHMM_trans_mat, try_hyperhMm$HyperHMM_trans_mat)
    expect_equal(try_correct$HyperHMM_trans_mat, try_hyperhmM$HyperHMM_trans_mat)
    expect_equal(try_correct$HyperHMM_trans_mat, try_HyperHmm$HyperHMM_trans_mat)
    expect_equal(try_correct$HyperHMM_trans_mat, try_HyperHMm$HyperHMM_trans_mat)
    expect_equal(try_correct$HyperHMM_trans_mat, try_HyperhMM$HyperHMM_trans_mat)
    expect_equal(try_correct$HyperHMM_trans_mat, try_HyperHmM$HyperHMM_trans_mat)
    expect_equal(try_correct$HyperHMM_trans_mat, try_hyperHMM$HyperHMM_trans_mat)
})


cat("\n Done test.HyperHMM.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
