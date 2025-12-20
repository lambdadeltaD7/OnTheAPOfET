library(knitr)
library(kableExtra)

get_table_caption <- function(filename) {
    name <- sub("\\.csv$", "", filename)
    parts <- strsplit(name, "_")[[1]]
    
    tests <- c(ET = "Energy", AD = "Anderson-Darling", CM = "Cramér-von Mises", 
               KS = "Kolmogorov-Smirnov", WMW = "Wilcoxon-Mann-Whitney")
    dists <- c(NORMAL = "Normal", CAUCHY = "Cauchy")
    
    test_full <- tests[parts[1]]
    dist_full <- dists[parts[2]]
    fixed_param <- parts[3]
    
    if(grepl("h1=", fixed_param)) {
        varying <- "h2"
    } else {
        varying <- "h1"
    }
    
    return(paste(test_full, "Test for", dist_full, "Distribution ( fixed", fixed_param, ")"))
}

for( fname in list.files(path = "./notes/res16/csvs", full.names = F)){
  df <- read.csv(paste("./notes/res16/csvs/", fname, sep=""))
  if(grepl("h1", fname)) {
        new_names <- gsub("h2\\.", "h2=", names(df))
    } else {
        new_names <- gsub("h1\\.", "h1=", names(df))
    }
  new_names[1] <- ""
  names(df) <- new_names
  table <- kable(df, "latex", booktabs = TRUE,
    caption = get_table_caption(fname), digits = 3) %>%
    kable_styling(latex_options = "hold_position")
  table <- gsub("\\\\begin\\{tabular\\}", "\\\\hspace*{-1.0cm}\\\\begin\\{tabular\\}\n", table)
  table <- gsub("\\\\addlinespace", "", table)
  table <- gsub("\\\\\\n\\AP", "\\\\\n\\\\addlinespace\nAP", table)
  writeLines(table, paste("./notes/res16/tables/table_et/", substr(fname ,1, nchar(fname)-4), ".tex" ,sep=""))
}