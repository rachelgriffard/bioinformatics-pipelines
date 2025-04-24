# iteratively create yaml files for tapestri
# Rachel Griffard-Smith
# 042425
# Last updated: 042425

library(yaml)

temp = yaml.load_file('config_template.yaml')

# Function to replace "sample" with sample names from a text file
replace_sample_names <- function(yaml_file, sample_names_file, output_file) {
  # Read YAML file
  yaml_content <- yaml.load_file(yaml_file)
  
  # Read sample names from text file
  sample_names <- readLines(sample_names_file)
  
  # Replace each instance of "sample" with corresponding sample name
  replace_sample_recursive <- function(obj, sample_name) {
    if (is.character(obj)) {
      return(gsub("sample", sample_name, obj))
    } else if (is.list(obj)) {
      return(lapply(obj, replace_sample_recursive, sample_name))
    } else {
      return(obj)
    }
  }
  
  # Iterate over sample names
  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    yaml_content <- replace_sample_recursive(yaml_content, sample_name)
  }
  
  # Write the modified YAML content to a new file
  write(yaml::as.yaml(yaml_content), output_file)
}

# Example usage
# Provide the paths to your YAML file, sample names file, and the output file
yaml_file <- "input.yml"
sample_names_file <- "sample_names.txt"
output_file <- "sample_config.yml"

replace_sample_names(yaml_file, sample_names_file, output_file)

cat("Replacement complete. Saved to", output_file, "\n")