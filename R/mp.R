mp <- function(mat, ...)
{
    image(max(mat, na.rm = TRUE) - t(mat[dim(mat)[1]:1, ]), ...)
}

