## assign each of the positions to a grid
## input is numbers, e.g. 3, 5, 10, 15
## and a windowSize, like 5
## output is 1-based on grid coordinates, like
## 1-5 -> 0, 6-10 -> 1, etc
## remove holes, sigma will be made able to handle it with bounding
assign_positions_to_grid <- function(
    L,
    gridWindowSize,
    method = "normal",
    smoothed_rate_cM = NULL,
    desired_gridWindowSize = NULL
) {
    if (method == "normal") {
        if (is.na(gridWindowSize) == FALSE) {
            grid <- ceiling(L / gridWindowSize)
            grid <- grid - min(grid)
            ## for L_grid, get first mid-point
            L_grid_start <- gridWindowSize * (ceiling(L[1] / gridWindowSize) - 0.5)
            grid_distances <- diff(unique(grid)) * gridWindowSize
            L_grid <- L_grid_start + c(0, cumsum(grid_distances))
            grid <- match(grid, unique(grid)) - 1
            nGrids <- length(grid_distances) + 1
            if ((length(L_grid) - length(grid_distances)) != 1) {
                stop("An error has been made assigning SNP positions to grid. Please report this")
            }
            ## this allows one to go from grid to start and end of SNPs in that grid
            ## so e.g. the fifth grid with 0-based index 4
            ## has 1-based SNPs starting from
            ## snps_in_grid_1_based[4 + 1, "grid_starts"]
            ## to
            ## snps_in_grid_1_based[4 + 1, "grid_ends"]        
            snps_in_grid_1_based <- cbind(
                snps_start = match(unique(grid), grid),
                snps_end = length(grid) - match(unique(grid), grid[length(grid):1]) + 1
            )
        } else {
            grid <- 0:(length(L) - 1)
            grid_distances <- diff(L)
            L_grid <- L
            nGrids <- length(L)
            snps_in_grid_1_based <- cbind(
                snps_start = 1:nGrids,
                snps_end = 1:nGrids
            )
        }
        return(
            list(
                grid = grid,
                grid_distances = grid_distances,
                L_grid = L_grid,
                nGrids = nGrids,
                snps_in_grid_1_based = snps_in_grid_1_based
            )
        )        
    } else if (method == "approximative") {
        grid <- attempt_to_better_grid(
            smoothed_rate_cM = smoothed_rate_cM,
            desired_gridWindowSize = desired_gridWindowSize,
            L = L
        )
        return(assign_grid_variables(grid, L))
    }
    return(to_out)
}
