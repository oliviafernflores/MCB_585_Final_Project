# Load libraries
library(tidyverse)
library(skimr)
library(rlang)
library(sva)

# Read raw data
file_path <- "/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/GSE13576_series_matrix.txt"
data <- readLines(file_path)

# The data is stored in a matrix, so I need to find the start of the line where
# the data actually starts after all the header information.
start_index <- which(str_detect(data, "^!series_matrix_table_begin"))[1]
if (!is.na(start_index)) {
  matrix_lines <- data[(start_index + 1):length(data)]
  # Same to find the ending of the matrix
  end_index <- which(str_detect(matrix_lines, "^!series_matrix_table_end"))[1]
  if (!is.na(end_index)) {
    matrix_lines <- matrix_lines[1:(end_index - 1)]
    # Now turn the matrix into a dataframe
    matrix_data <- str_split(matrix_lines, "\t") %>% 
      map(~str_trim(.)) %>%
      map(~na_if(.,""))
    df <- as.data.frame(matrix_data, stringsAsFactors = FALSE)
    colnames(df) <- df[1,]
    df <- df[-1,]
  } else {
    print("No matrix end")
  }
} else {
  print("No matrix start")
}

# Debugging to check that the data was actually turned into a dataframe
df <- as.data.frame(matrix_data, stringsAsFactors = FALSE)

# Fix the column names
colnames(df) <- df[1,]
df <- df[-1, ]

# Fix the index
df <- df %>%
  column_to_rownames(var = '"ID_REF"') 


# Remove extra " from column names
df <- df %>%
  rename_with(~str_replace_all(., '"', ''))

# Remove extra " from row names
df$ID_REF <- str_replace_all(df$ID_REF, '"', '')

# Create a vector to store the sample types
sample_type <- c()

# Loop through the data to find the line with '!Sample_source_name_ch1' to get a vector of all the sample types
for (i in seq_along(data)) {
  if (grepl("^!Sample_source_name_ch1", data[i])) {
    sample_type <- str_split(data[i], "\t")[[1]]
    break
  }
}

# Remove the added label
sample_type <- sample_type[-1]

# Format the sample types in the sample_type vector
sample_type <- str_trim(sample_type)
sample_type <- gsub('"', '', sample_type)
sample_type <- gsub(" ", "_", sample_type)

# Add the sample types to the dataframe with the data
df$Sample_Type <- sample_type


# Do the same process for time status
time_status <- c()


for (i in seq_along(data)) {
  if (grepl("^!Sample_title", data[i])) {
    time_status <- str_split(data[i], "\t")[[1]]
    break
  }
}

time_status <- time_status[-1]

time_status <- sapply(time_status, function(time) {
  time <- str_trim(time) 
  time <- gsub('"', '', time)
  
  if (grepl("no relapse", time)) {
    return("none")
  } else if (grepl("early relapse", time)) {
    return("early")
  } else if (grepl("late relapse", time)) {
    return("late")
  } else if (grepl("xeno_loTTL", time)) {
    return("long_TTL")
  } else if (grepl("xeno_shTTL", time)) {
    return("short_TTL")
  } else {
    return(0)
  }
})

df$Time_Status <- time_status

# Now that the time info and sample type is added to the data, separate the data into the xenograft derived samples and the direct patient samples

df_xenograft <- df %>% filter(str_detect(Sample_Type, regex("xenograft", ignore_case = TRUE)))

df_patient <- df %>% filter(!str_detect(Sample_Type, regex("xenograft", ignore_case = TRUE)))

# A function to turn the matrix (all the numbers are currently stored as strings) into a matrix of numbers so that I can calculate differential gene expression
convert_strings_to_numeric <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      converted <- suppressWarnings(as.numeric(x))
      return(converted)
    } else {
      return(x)
    }
  })
  
  return(df)
}


# A function for the normalization
normalize_numeric_columns <- function(df, ref_df) {
  # save the rownames to be used if needed later
  original_rownames <- rownames(df)
  
  # Convert the numerical strings to numbers using the previous function
  df <- convert_strings_to_numeric(df) 
  
  # Separate the numeric columns
  numeric_cols <- df %>% select(where(is.numeric))
  non_numeric_cols <- df %>% select(where(~!is.numeric(.)))
  
  # Use scale() to normalize the numeric columns
  normalized_data <- scale(numeric_cols)
  
  # Make a new dataframe with the normalized data
  df_normalized <- cbind(non_numeric_cols, as.data.frame(normalized_data))
  
  # Restore the original rownames
  rownames(df_normalized) <- original_rownames
  
  # Restore column names
  colnames(df_normalized) <- c(colnames(non_numeric_cols), colnames(numeric_cols))
  
  return(df_normalized)
}



# Do normalization for the 'xenograft' and 'patient' dataframes
df_xenograft_normalized <- normalize_numeric_columns(df_xenograft, df)
df_xenograft_normalized$ID_REF <- c(df_xenograft$ID_REF, "GENE")
df_xenograft_normalized$Sample_Type <- c(sample_type[198:209], 'gene_info')
df_xenograft_normalized$Time_Status <- c(time_status[198:209], 'gene_info')
df_patient_normalized <- normalize_numeric_columns(df_patient, df)
df_patient_normalized$ID_REF <- c(df_patient$ID_REF, "GENE")
df_patient_normalized$Sample_Type <- c(sample_type[1:197], 'gene_info', 'gene_info')
df_patient_normalized$Time_Status <- c(time_status[1:197], 'gene_info', 'gene_info')

# Save the normalized data to CSV files
write_csv(df_xenograft_normalized, "/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/xenograft_data_normalized.csv")
write_csv(df_patient_normalized, "/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/patient_data_normalized.csv")
