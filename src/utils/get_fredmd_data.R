# Load libraries
library(BVAR)

# Load data
get_fredmd_data <- function(start_date) {
    # Load data
    data <- read.csv("data/current.csv", row.names = 1)
    transformations <- as.numeric(data[1, ])
    data <- data[-1, ]
    rownames(data) <- as.Date(rownames(data), format = "%m/%d/%Y")

    # Transform and clean data
    fred_md_transformed <- fred_transform(data, codes = transformations, na.rm = FALSE)
    fred_md_transformed <- fred_md_transformed[
        rownames(fred_md_transformed) >= start_date,
    ]
    fred_md_transformed <- fred_md_transformed[
        ,
        colSums(is.na(fred_md_transformed)) == 0
    ]
}
