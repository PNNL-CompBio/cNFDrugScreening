library("ggplot2")
library("pracma")
library("RGCCA")
library("cowplot")
library("ggpubr")

nan_norm <- function(x) {
  return(x / sqrt(sum(x[!is.na(x)]^2)))
}

phospho <- read.table(
  "data/merged_phospho.txt.gz",
  header = TRUE, 
  sep = ",",
  row.names = 1
)
drug_response <- read.table(
  "data/drug_responses.csv",
  sep = ",",
  header = TRUE, 
  row.names = 1
)

# Trim to patients with both phospho and drug response measurements
patients = intersect(rownames(phospho), rownames(drug_response))
phospho = phospho[patients,]
drug_response = drug_response[patients,]

blocks <- list(
  phospho = t(apply(phospho, 1, nan_norm)),
  drug_response = t(apply(drug_response, 1, nan_norm))
)
taus <- linspace(0, 1, 11)
sparsities <- linspace(0.25, 1, 7)

ranks = 6
plots <- list()
for (tau in taus) {
  for (sparsity in sparsities) {
    res <- rgcca(
      blocks, 
      ncomp = ranks,
      scheme = "factorial",
      method = "sgcca",
      scale = TRUE,
      scale_block = TRUE,
      sparsity = sparsity,
      tau = tau,
      verbose = FALSE
    )
    
    aves <- data.frame(
      list(
        res$AVE$AVE_inner,
        res$AVE$AVE_outer,
        res$AVE$AVE_X$phospho,
        res$AVE$AVE_X$drug_response,
        1:ranks
      )
    )
    colnames(aves) = c(
      "Inner", "Outer", "Phospho", "Drug Response", "Components"
    )
    
    ave_plot <- ggplot(data = aves, aes(x = Components))
    ave_plot = ave_plot + geom_line(aes(y = Inner, color = "Inner"))
    ave_plot = ave_plot + geom_line(aes(y = Outer, color = "Outer"))
    ave_plot = ave_plot + geom_line(aes(y = Phospho, color = "Phospho"))
    ave_plot = ave_plot + geom_line(
      aes(y = `Drug Response`, color = "Drug Response")
    )
    ave_plot = ave_plot + ylab("AVE")
    ggsave(
      sprintf(
        "output/tau_%d_sparsity_%d.png",
        as.integer(tau * 100), 
        as.integer(sparsity * 100)
      )
    )
  }
}

res <- rgcca(
  blocks, 
  ncomp = 4,
  scheme = "factorial",
  method = "sgcca",
  scale = TRUE,
  scale_block = TRUE,
  sparsity = sparsity,
  tau = tau,
  verbose = FALSE
)
write.csv(res$a$phospho, "output/phospho_loadings.csv")
write.csv(res$Y$phospho, "output/phospho_scores.csv")
write.csv(res$a$drug_response, "output/drug_loadings.csv")
write.csv(res$Y$drug_response, "output/drug_scores.csv")
