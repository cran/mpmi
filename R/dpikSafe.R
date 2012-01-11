# Function to handle cases where dpik() is unable to estimate
# a bandwidth, usually because a data vector has an interquartile range of 0.
dpikSafe <- function(x, ...)
{
	result <- try(dpik(x, ...), silent = TRUE)

	if (class(result) == "try-error")
	{
		msg <- geterrmessage()
		if (grepl("scale estimate is zero for input data", msg))
		{
			warning("Using standard deviation as scale estimate, probably because IQR == 0")
			result <- dpik(x, scalest = "stdev", ...)	
		} else 
		{
			stop(msg)
		}		
		
	}
	return(result)
}
