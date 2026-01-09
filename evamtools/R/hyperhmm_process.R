translate_state <- function(n_decimal, genes) { 
  L <- length(genes) 
  if (n_decimal == 0) return("WT") 
  bits <- as.integer(intToBits(n_decimal)) # Convertir decimal a posiciones de bits activos 
  positions <- which(bits[1:L] == 1)   # Encontrar qué posiciones tienen un 1 (ajustando al orden de HyperHMM) 
  # Mapear posiciones a nombres de genes 
  mutated_names <- sort(genes[L - positions + 1])   # Nota: HyperHMM cuenta de izquierda a derecha, R de derecha a izquierda 
  return(paste(mutated_names, collapse = ", ")) 
} 

do_HyperHMM <- function(data,
                        opts = list(initialstates = NULL, 
                                    seed = 1L, 
                                    nboot = 100L, 
                                    fullsample = 1L, 
                                    outputinput = 0L),
                        silent = FALSE) {

    feature_labels <- colnames(data)
    num_features <- ncol(data)
    if(!is.matrix(data)) data <- as.matrix(data)

    #-------- AJUSTAMOS EL MODELO
    args <- c(list(obs = data), opts)
    out <- do.call(hyperhmm::HyperHMM, args)
    # Filtramos por bootstrap = 0
    boot0 <- out$transitions[out$transitions$Bootstrap == 0, ]

    #--------- IDs GENOTIPOS
    all_ids <- 0:(2^(num_features)-1)
    decoded_states <- vapply(all_ids, translate_state,
                            character(1),
                            genes = feature_labels)
    num_muts <- vapply(strsplit(decoded_states, ", "), function(x) {
        if (length(x) == 1 && x == "WT") return(0)
        return(length(x))
    }, numeric(1))
    
    ordered_names <- decoded_states[order(num_muts, decoded_states)]
    genotype_id_ordered <- setNames(1:length(ordered_names), ordered_names)

    #--------- TABLA EDGES
    from_names <- vapply(boot0$From, translate_state, character(1), genes = feature_labels)
    to_names <- vapply(boot0$To, translate_state, character(1), genes = feature_labels)
    edges <- data.frame(From = from_names,
                        To = to_names,
                        Probability = boot0$Probability,
                        Relation = "Single",
                        stringsAsFactors = FALSE)
    edges$Edges = paste (edges$From,edges$To, sep = " -> ")

    #-------- TRANSITION MATRIX
    trans_mat <- Matrix(0,
                        nrow = length(ordered_names),
                        ncol = length(ordered_names),
                        sparse = TRUE)
    rownames(trans_mat) <- colnames(trans_mat) <- ordered_names

    row_idx <- match(edges$From, ordered_names)
    col_idx <- match(edges$To, ordered_names)
    
    trans_mat[cbind(row_idx, col_idx)] <- edges$Probability
    trans_mat <- as(trans_mat, "generalMatrix") #Forzamos conversión a dgCMatrix

    #-------- TRANS FLUX MAT
    trans_rate_mat <- Matrix(0,
                        nrow = length(ordered_names),
                        ncol = length(ordered_names),
                        sparse = TRUE)
    rownames(trans_rate_mat) <- colnames(trans_rate_mat) <- ordered_names

    row_idx_flux <- match(edges$From, ordered_names)
    col_idx_flux <- match(edges$To, ordered_names)

    trans_rate_mat[cbind(row_idx_flux, col_idx_flux)] <- boot0$Flux
    trans_rate_mat <- as(trans_rate_mat, "dgCMatrix") #Forzamos conversión a dgCMatrix


    #--------- FRECUENCIAS PREDICHAS DE LOS GENOTIPOS
    # uses the flux (effective flux) to get a more accurate representation 
    # of the system's evolutionary eq and calculate the predicted genotype freqs
    pre_genotype_freqs <- probs_from_trm(trans_rate_mat, all_genotypes = TRUE)
    
    features=list(N=out$L, names_features=feature_labels)

    return(list(
        trans_matrix = trans_mat,
        predicted_freqs = pre_genotype_freqs,
        trans_rate_mat = trans_rate_mat,
        stats = out$stats,
        paths_all = out$viz,
        n_features = features,
        raw_output = out
    ))
}