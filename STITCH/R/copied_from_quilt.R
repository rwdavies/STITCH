modified_calculate_pse <- function(
    test,
    truth,
    LL,
    seed = NULL
) {
    ## for testing
    if (is.null(seed) == FALSE)
        set.seed(seed)
    rownames(test) <- 1:nrow(test)
    which_sites <-
        rowSums(truth == 0 | truth == 1, na.rm = TRUE) == 2 &
        rowSums(truth, na.rm = TRUE) == 1 &
        rowSums(is.na(truth)) == 0
    truth <- truth[which_sites, ]
    test <- test[which_sites, ]
    if (nrow(test) == 0)
        return(NA)
    ## as these sites, discrepency as well
    disc <- sum(round(rowSums(test)) != 1)
    test <- round(test)
    ## round test data for now. choose hets at random
    ## for homs
    test[, 1] <- as.integer(round(test[, 1]))
    test[, 2] <- as.integer(round(test[, 2]))
    testO <- test
    truthO <- truth
    for(i_option in 1:2) {
        test <- testO
        truth <- truthO
        ## specifically remove from consideration double phase switch errors
        w <- rowSums(test) == 1
        w2 <- which(diff(abs(test[w,1] - truth[w,1])) != 0)
        to_remove <- NULL
        if (length(w2) > 0) {
            w3 <- which(diff(w2) == 1)
            if (length(w3) > 0) {
                for(a in w3) {
                    c <- w2[c(a, a + 1)]
                    to_remove <- c(to_remove, which(w)[c])
                }
            }
        }
        ## double pse are two consecutive
        if (length(to_remove) > 0) {
            test <- test[-to_remove, ]
            truth <- truth[-to_remove, ]
        }
        ##
        if (i_option == 1) {
            ## only consider non-discrepent sites
            ## chose best start
            if (test[1, 1] != truth[1, 1])
                test <- test[, c(2, 1)]
            ## calculate number of differences
            w <- rowSums(test) == 1
            if (sum(w) == 0) {
                print("Test has no hets! possibly an error or homo over region, possibly no record dosages turned on in impute_all")
                switches1 <- cbind(i1 = NA, i2 = NA, l1 = NA, l2 = NA)
                phase_errors_def1 <- 0
                phase_sites_def1 <- 0
            } else {
                y <- diff(abs(test[w,1] - truth[w,1])) != 0
                phase_errors_def1 <- sum(y)
                phase_sites_def1 <- sum(w) - 1
                s <- as.integer(rownames(test[w, , drop = FALSE][c(as.logical(y), FALSE), , drop = FALSE]))
                e <- as.integer(rownames(test[w, , drop = FALSE][c(FALSE, as.logical(y)), , drop = FALSE]))
                switches1 <- cbind(i1 = s, i2 = e, l1 = LL[s], l2 = LL[e])
            }
        }
        if (i_option == 2) {
            choose_at_random <- which(rowSums(test) != 1)
            if (length(choose_at_random) > 0) {
                test[choose_at_random, ] <- 0
                r <- sample(
                    c(1, 2),
                    length(choose_at_random),
                    replace = TRUE
                )
                test[cbind(choose_at_random, r)] <- 1
            }
            ## chose best start
            if (test[1, 1] != truth[1, 1])
                test <- test[, c(2, 1)]
            ## calculate number of differences
            phase_errors_def2 <- sum(diff(abs(test[,1] - truth[,1])) != 0)
            phase_sites_def2 <- nrow(test) - 1
        }
    }
    ##
    return(
        list(
            values = c(
                phase_errors_def1 = phase_errors_def1,
                phase_sites_def1 = phase_sites_def1,
                phase_errors_def2 = phase_errors_def2,
                phase_sites_def2 = phase_sites_def2,
                disc_errors = disc,
                dist_n = nrow(test)
            ),
            switches1 = switches1
        )
    )
}
