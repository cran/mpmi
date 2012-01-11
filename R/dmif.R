dmi.pw <- function(disc1, disc2)
{
    disc1 <- as.integer(factor(disc1))
    disc2 <- as.integer(factor(disc2))

    ok <- is.finite(disc1) & is.finite(disc2)

    disc1 <- disc1[ok]
    disc2 <- disc2[ok]

    mi <- 0.0
    bcmi <- 0.0
    zvalue <- 0.0

    result <- .Fortran("dmi",
                       disc1,
                       as.integer(length(disc1)),
                       disc2,
                       as.integer(length(disc2)),
                       mi = mi, 
                       bcmi = bcmi,
                       zvalue = zvalue,
                       DUP = TRUE)

    return(result[c("mi", "bcmi", "zvalue")])
}

dmi <- function(dmat)
{
    # Convert to integers
    dint <- matrix(0L, nrow = nrow(dmat), ncol = ncol(dmat))
    for (i in 1: ncol(dmat))
    {
        dint[,i] <- as.integer(factor(dmat[,i]))
    }

    mi <- matrix(0.0, nrow = ncol(dint), ncol = ncol(dint))
    bcmi <- matrix(0.0, nrow = ncol(dint), ncol = ncol(dint))
    zans <- matrix(0.0, nrow = ncol(dint), ncol = ncol(dint))

    result <- .Fortran("dmim",
                       dint,
                       as.integer(nrow(dint)),
                       as.integer(ncol(dint)),
                       mi = mi, 
                       bcmi = bcmi, 
                       zvalues = zans,
                       NAOK = TRUE, 
                       DUP = TRUE)

    return(result[c("mi", "bcmi", "zvalues")])
}

dminjk.pw <- function(disc1, disc2)
{
    disc1 <- as.integer(factor(disc1))
    disc2 <- as.integer(factor(disc2))
    ans <- as.double(0)
    return(.Fortran("dminjk",
                    disc1,
                    as.integer(length(disc1)),
                    disc2,
                    as.integer(length(disc2)),
                    result = ans, DUP = TRUE)$result)
}

dminjk <- function(dmat)
{
    # Convert to integers
    dint <- matrix(0L, nrow = nrow(dmat), ncol = ncol(dmat))
    for (i in 1: ncol(dmat))
    {
        dint[,i] <- as.integer(factor(dmat[,i]))
    }

    ans <- matrix(0.0, nrow = ncol(dint), ncol = ncol(dint))
    return(.Fortran("dmimnjk",
                    dint,
                    as.integer(nrow(dint)),
                    as.integer(ncol(dint)),
                    result = ans, NAOK = TRUE, DUP = TRUE)$result)
}

