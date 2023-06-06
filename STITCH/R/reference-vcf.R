## a simple hacky way to do this, to check Rcpp functionality mostly
R_get_hap_info_from_vcf <- function(
    vcffile,
    af_cutoff,
    region
) {
    ## write data parser
    ## input: phased vcf, maf threshold
    ## output:
    ##  - at common sites: rhb_t, then later maybe, zilong indices, robbie things (e.g. hapMatcher)
    ##  - at rare sites: an object that somehow can be used to rebuild
    cmd <- paste0("bcftools view ", shQuote(vcffile), " ", region)
    haps <- data.table::fread(cmd = cmd, stringsAsFactors = FALSE, data.table = FALSE)
    ## skip if don't meet the conditions
    keep <-
        (nchar(haps[, 4]) == 1) &
        (nchar(haps[, 5]) == 1)
    ## remove second (and more) if position the same
    keep[ which(diff(haps[, 2]) == 0) + 1] <- FALSE
    n_skipped <- sum(!keep)
    haps <- haps[keep, ]
    ## now do stuff    
    pos <- haps[, c(2, 4, 5)] ## simpler, add back in later
    rownames(pos) <- 1:nrow(pos)
    LAll <- haps[, 2]
    h <- haps[, -(1:9)]
    h <- as.matrix(h)
    K <- ncol(h) * 2
    nSNPsAll <- nrow(h)
    rhi <- array(as.integer(NA), c(nSNPsAll, K)) ## reference haplotype integers
    rhi[, seq(1, K, 2)] <- matrix(as.integer(substr(h, 1, 1)), nrow = nSNPsAll)
    rhi[, seq(2, K, 2)] <- matrix(as.integer(substr(h, 3, 3)), nrow = nSNPsAll)
    x <- rowSums(rhi)
    ref_alleleCount <- cbind(x, K, x / K)
    colnames(ref_alleleCount) <- NULL
    ## filtering 
    snp_is_common <- !(ref_alleleCount[, 3] < af_cutoff)
    L <- LAll[snp_is_common]
    nSNPs <- sum(snp_is_common)
    nGrids <- ceiling(nSNPs / 32)
    grid <- as.integer(floor((1:nSNPs) / 32))
    snp_is_common_1_based <- which(snp_is_common)
    ## not sure of the right way to do this on the fly? with rhb or rhb_t?
    rhb <- matrix(0L, nrow = nGrids, ncol = K)
    for(iGrid in 1:nGrids) {
        s <- min(32 * (iGrid - 1) + 1, nSNPs)
        e <- min(32 * iGrid, nSNPs)
        w <- snp_is_common_1_based[s:e]
        for(k in 1:K) {
            rhb[iGrid, k] <- rcpp_int_contract(rhi[w, k])
        }
    }
    rhb_t <- t(rhb)
    ##
    ##
    snp_is_rare_1_based <- which(!snp_is_common)    
    rare_per_hap_info <- sapply(1:K, function(k) {
        ## for now store position among all SNPs
        snp_is_rare_1_based[which(rhi[!snp_is_common, k] == 1)]
    })
    ##
    return(
        list(
            pos = pos,
            rhb_t = rhb_t,
            rare_per_hap_info = rare_per_hap_info,
            snp_is_common = snp_is_common,
            ref_alleleCount = ref_alleleCount,
            n_skipped = n_skipped
        )
    )
}
