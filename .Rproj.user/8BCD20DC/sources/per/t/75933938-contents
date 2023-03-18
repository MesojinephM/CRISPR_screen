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
                as.factor
                )
         )

# View data summary
summary(mageck)

# Sort genes by `lfc` and create named vector of `lfc` values
lfc_sort <- mageck[order(mageck$lfc, decreasing = TRUE),]
mageck_lfc_sort <- as.vector(lfc_sort$lfc)
names(mageck_lfc_sort) <- lfc_sort$id
head(mageck_lfc_sort)


######################################################################################

# Import gene set lists
# Code adapted from https://stackoverflow.com/questions/6602881/text-file-to-list-in-r

# Read in pathways data as list and split elements into strings
gmt <- scan("m2.cp.v2022.1.Mm.symbols.gmt", what = "", sep = "\n")
pathways <- strsplit(gmt, "[[:space:]]+")

# Assign first entry to names of each list element
names(pathways) <- sapply(pathways, `[[`, 1)

# save urls to separate list for reference
source <- sapply(pathways, `[[`, 2)

# Change gene names to lowercase to reduce errors
pathways <- lapply(pathways[-1], tolower)

# Remove two beginning reference rows
pathways <- lapply(pathways, `[`, -c(1:2))

# Preview data set
head(pathways)


###############################################################################

# Run GSEA on pre-ranked gene expression data
score_df <- fgsea::fgseaSimple(
  pathways = pathways,
  stats = mageck_lfc_sort,
  nperm = 1000,
  minSize = 5,
  scoreType = "std",
  gseaParam = 1
)

# Add variable containing number of genes in the leading edge subset
# score_df$nLeadingEdge <- as.numeric(lapply(score_df$leadingEdge, length))

# Write GSEA results to tab-delimited file, excluding leading edge variable
write.table(
  x = score_df[,1:7,9],
  file = "gsea_output.txt",
  sep = "\t", 
  row.names = FALSE
)

# Create leading edge list object and write to external file
leading_edge <- score_df$leadingEdge
names(leading_edge) <- score_df$pathway

utils::capture.output(leading_edge, 
                      file = "leading_edge.csv") 
utils::capture.output(leading_edge, 
                      file = "leading_edge.txt") 

# Store GSEA output and data objects as .RData file
save(score_df,
     mageck_lfc_sort,
     pathways,
     file = "gsea_output.RData"
     )

# FGSEA citation
# citEntry(entry="ARTICLE",
#          author = c(person("Gennady", "Korotkevich"),
#                     person("Vladimir", "Sukhov"),
#                     person("Alexey", "Sergushichev")),
#          title = "Fast gene set enrichment analysis",
#          year = "2019",
#          doi = "10.1101/060012",
#          publisher = "Cold Spring Harbor Labs Journals",
#          URL = "http://biorxiv.org/content/early/2016/06/20/060012",
#          journal = "bioRxiv",
#          textVersion = paste("G. Korotkevich, V. Sukhov, A. Sergushichev.",
#                              "Fast gene set enrichment analysis.",
#                              "bioRxiv (2019), doi:10.1101/060012")
# )