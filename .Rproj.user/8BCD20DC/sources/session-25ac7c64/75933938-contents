#########################
# GSEA Code             #
# Marie Moriarty        #
# Senior Capstone       #
# Spring 23             #
#########################

# Import MAGeCK data set
library(readr)
library(ggplot2)
library(dplyr)

data <- read_delim("mageckRRA.gene_summary.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Reassign column names
colnames(data) <- c("id",
                      "num",
                      "neg.score",
                      "neg.p_value",
                      "neg.fdr",
                      "neg.rank",
                      "neg.goodsgrna",
                      "neg.lfc",
                      "pos.score",
                      "pos.p_value",
                      "pos.fdr",
                      "pos.rank",
                      "pos.goodsgrna",
                      "lfc"
)

# Filter and clean data
mageck <- data %>%
  # Remove observations where `num` for sgRNA != 4
  filter(num == "4") %>%
  # Drop unnecessary `num` and `neg.lfc` columns.
  select(-c(neg.lfc, num)) %>%
  # Change `id` column to all lowercase characters.  
  mutate(across(id, tolower)) %>%
  # Convert to type factor
  mutate(across(c(neg.goodsgrna,
                  pos.goodsgrna),
                as.factor))


# view data summary
summary(mageck)
