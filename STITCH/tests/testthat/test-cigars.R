## base qualities from 25 to 35 inclusive
## paste0(sapply(33 + 25 + 0:9, function(x) rawToChar(as.raw(x))), collapse = "")
## :;<=>?@ABC

test_that("can calculate read length with one element", {
    readLength <- get_read_span(1, "M")
    expect_equal(readLength, 1)
})

test_that("can calculate read length with two elements", {
    readLength <- get_read_span(c(1, 2), c("M", "M"))
    expect_equal(readLength, 3)
})

test_that("can calculate read length with ten elements", {
    readLength <- get_read_span(rep(3, 10), rep("M", 10))
    expect_equal(readLength, 30)
})

test_that("can calculate read length using M and D", {
    readLength <- get_read_span(c(5, 6), c("M", "D"))
    expect_equal(readLength, 11)
})

test_that("can calculate read length skipping I", {
    readLength <- get_read_span(c(5, 6, 7, 8), c("M", "D", "I", "M"))
    expect_equal(readLength, 5 + 6 + 8)
})


test_that("can split cigar strings", {
    cigarRead <- c("101M", "101M")
    splitCigarRead <- cpp_cigar_split_many(cigarRead)
    expect_equal(    
        splitCigarRead,
        list(
            list(101, "M"),
            list(101, "M")
        )
    )
})

test_that("complicated cigar strings are split properly", {
    cigarRead <- c("11M12I13M", "14M15D16M")
    splitCigarRead <- cpp_cigar_split_many(cigarRead)
    expected_result <- list(
        list(c(11, 12, 13), c("M", "I", "M")),
        list(c(14, 15, 16), c("M", "D", "M"))
    )
    expect_equal(splitCigarRead, expected_result)
})

test_that("a mixture of complicated and uncomplicated cigar strings are split properly", {
    cigarRead <- c("11M12I13M", "14M15D16M", "101M", "102M")
    splitCigarRead <- cpp_cigar_split_many(cigarRead)
    expected_result <- list(
        list(c(11, 12, 13), c("M", "I", "M")),
        list(c(14, 15, 16), c("M", "D", "M")),
        list(101, "M"),
        list(102, "M")
    )
    expect_equal(splitCigarRead, expected_result)
})


## soft clipped read
## "pos" given is first not-soft clipped base
test_that("a soft clipped read is filtered appropriately when not using soft clipped bases", {
    splitCigarRead <- list(list(c(3, 3, 4), c("S", "M", "S")))
    out <- deal_with_soft_clipped_bases(
        splitCigarRead = splitCigarRead,
        useSoftClippedBases = FALSE,
        posRead = 10,
        seqRead = "AAACCCGGGG",
        qualRead = ":;<=>?@ABC"
    )
    expected_result <- list(list(3, "M"))
    expect_equal(out$splitCigarRead, expected_result)
    expect_equal(out$seqRead, "CCCGGGG")
    expect_equal(out$qualRead, "=>?@ABC")
    expect_equal(out$posRead, 10)
})

test_that("a soft clipped read is filtered appropriately when with a normal read when not using soft clipped bases", {
    splitCigarRead <- list(
        list(c(21, 22, 23), c("M", "I", "M")),
        list(c(1, 1, 1), c("S", "M", "S")),
        list(c(14, 15, 16), c("M", "D", "M"))
    )
    out <- deal_with_soft_clipped_bases(
        splitCigarRead = splitCigarRead,
        useSoftClippedBases = FALSE,
        posRead = c(5, 10, 15),
        seqRead = c("not_used", "ACG", "not_used"),
        qualRead = c("not_used", ":;<", "not_used")
    )
    expected_result <- list(
        list(c(21, 22, 23), c("M", "I", "M")),
        list(1, "M"),
        list(c(14, 15, 16), c("M", "D", "M"))
    )
    expect_equal(out$splitCigarRead, expected_result)
    expect_equal(out$seqRead, c("not_used", "CG", "not_used"))
    expect_equal(out$qualRead, c("not_used", ";<", "not_used"))
    expect_equal(out$posRead, c(5, 10, 15))
})

test_that("a soft clipped read is filtered appropriately when using soft clipped bases", {
    splitCigarRead <- list(list(c(2, 3, 4), c("S", "M", "S")))
    out <- deal_with_soft_clipped_bases(
        splitCigarRead = splitCigarRead,
        useSoftClippedBases = TRUE,
        posRead = 3,
        seqRead = "AACCCGGGG",
        qualRead = ":;<=>?@ABC"
    )
    expected_result <- list(list(c(2, 3, 4), c("M", "M", "M")))
    expect_equal(out$splitCigarRead, expected_result)
    expect_equal(out$seqRead, "AACCCGGGG")
    expect_equal(out$qualRead, ":;<=>?@ABC")
    expect_equal(out$posRead, 3 - 2)
})

test_that("a soft clipped read is filtered appropriately when with a normal read when using soft clipped bases", {
    splitCigarRead <- list(
        list(c(11, 12, 13), c("M", "I", "M")),
        list(c(3, 2, 1), c("S", "M", "S")),
        list(c(14, 15, 16), c("M", "D", "M"))
    )
    out <- deal_with_soft_clipped_bases(
        splitCigarRead = splitCigarRead,
        useSoftClippedBases = TRUE,
        posRead = c(5, 10, 15),
        seqRead = c("not_used", "AAACCG", "not_used"),
        qualRead = c("not_used", ":;<=>?", "not_used")
    )
    expected_result <- list(
        list(c(11, 12, 13), c("M", "I", "M")),
        list(c(3, 2, 1), c("M", "M", "M")),
        list(c(14, 15, 16), c("M", "D", "M"))
    )
    expect_equal(out$splitCigarRead, expected_result)
    expect_equal(out$seqRead, c("not_used", "AAACCG", "not_used"))
    expect_equal(out$qualRead, c("not_used", ":;<=>?", "not_used"))
    expect_equal(out$posRead, c(5, 10 - 3, 15))
})
