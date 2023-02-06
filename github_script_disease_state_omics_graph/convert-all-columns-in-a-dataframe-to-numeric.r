df[] <- lapply(df, function(x) as.numeric(as.character(x)))
df
sapply(df, class)