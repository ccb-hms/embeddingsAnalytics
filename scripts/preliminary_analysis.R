# Load required packages
library(SingleCellExperiment)
library(scater)
library(scDiagnostics)

# ___________________
# Data Preprocessing
# ___________________

# Read in data
full_data <- read.csv("data/bge-small-en-v1.5_embedding.csv")

# Format corpus embeddings
corpus_embeddings <- matrix(NA, nrow = nrow(full_data), ncol = 384)
file_name <- vector(length = nrow(full_data))
page_number <- numeric(length = nrow(full_data))

for(i in 1:nrow(full_data)){
  
  char_string <- gsub("\\[|\\]", "", full_data$Node.Embedding[i])  
  num_strings <- strsplit(char_string, ", ")[[1]] 
  corpus_embeddings[i,] <- as.numeric(num_strings)
  file_name[i] <- full_data$Filename[i]
  page_number[i] <- full_data$Page.No.[i]
}

colnames(corpus_embeddings) <- paste0("X", 1:ncol(corpus_embeddings))
corpus_embeddings <- data.frame(data_type = rep("corpus", nrow(corpus_embeddings)), 
                                file_name = file_name, 
                                page_number = page_number, 
                                corpus_embeddings)
corpus_embeddings <- corpus_embeddings[seq(1, nrow(corpus_embeddings), by = 2),]
rownames(corpus_embeddings) <- paste0("Corpus", 1:nrow(corpus_embeddings))

# Format question embeddings
question_embeddings <- matrix(NA, nrow = nrow(full_data), ncol = 384)
file_name <- vector(length = nrow(full_data))
page_number <- numeric(length = nrow(full_data))

for(i in 1:nrow(full_data)){
  
  char_string <- gsub("\\[|\\]", "", full_data$Question.Embedding[i])  
  num_strings <- strsplit(char_string, ", ")[[1]] 
  question_embeddings[i,] <- as.numeric(num_strings)
  file_name[i] <- full_data$Filename[i]
  page_number[i] <- full_data$Page.No.[i]
}

colnames(question_embeddings) <- paste0("X", 1:ncol(question_embeddings))
rownames(question_embeddings) <- paste0("Question", 1:nrow(question_embeddings))
question_embeddings <- data.frame(data_type = rep("question", nrow(question_embeddings)), 
                                  file_name = file_name, 
                                  page_number = page_number, 
                                  question_embeddings)

# ___________________________________
# Multivariate Analysis (File Level)
# ___________________________________

# Create SCE objects from embeddings
corpus_experiment <- SingleCellExperiment(
  assays = list(logcounts = t(corpus_embeddings[, paste0("X", 1:384)])), 
  colData = corpus_embeddings[, c("file_name", "page_number")])
questions_experiment <- SingleCellExperiment(
  assays = list(logcounts = t(question_embeddings[, paste0("X", 1:384)])), 
  colData = question_embeddings[, c("file_name", "page_number")])

# File types
file_names <- c("AMIE.pdf", 
                "ImageNet.pdf", "Unet.pdf", # Two CNN papers with similar topics
                "PRISM.pdf", 
                "missForestMICE.pdf")

# Subset embedding experiments
corpus_experiment <- corpus_experiment[, which(corpus_experiment$file_name %in% file_names)] 
questions_experiment <- questions_experiment[, which(questions_experiment$file_name %in% file_names)] 

# Generate the MDS scatter plot with cell type coloring
mds_plot <- plotCellTypeMDS(query_data = questions_experiment, 
                            reference_data = corpus_experiment, 
                            cell_types = file_names,
                            query_cell_type_col = "file_name", 
                            ref_cell_type_col = "file_name")
print(mds_plot)

# Plot the PC data
corpus_experiment <- runPCA(corpus_experiment)
questions_experiment <- runPCA(questions_experiment)
pc_plot <- plotCellTypePCA(query_data = questions_experiment, 
                           reference_data = corpus_experiment,
                           cell_types = file_names,
                           query_cell_type_col = "file_name", 
                           ref_cell_type_col = "file_name", 
                           pc_subset = c(1:10))
print(pc_plot)

# Plot discriminant space
disc_output <- calculateDiscriminantSpace(reference_data = corpus_experiment,
                                          query_data = questions_experiment, 
                                          query_cell_type_col = "file_name", 
                                          ref_cell_type_col = "file_name",
                                          eigen_threshold  = 1e-15,
                                          n_tree = 500,
                                          n_top = 20,
                                          calculate_metrics = FALSE,
                                          alpha = 0.01)
disc_plot <- plot(disc_output, plot_type = "scatterplot")
print(disc_plot)
