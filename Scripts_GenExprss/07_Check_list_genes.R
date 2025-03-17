# LOAD DATA USING THE 07_Enrichment_AvsB

# genes 
genes <- c("ATF4", "CHAC1", "SESN2", "ALDH1L2", "CBS", "PHGDH", "AMD1")
  
noDXvsDOX24H <- bg[which(bg$Comparison == "noDOXvsDOX24H"), ]
noDXvsDOX42H <- bg[which(bg$Comparison == "noDOXvsDOX42H"), ]
DOX24HvsDOX42H <- bg[which(bg$Comparison == "DOX24HvsDOX42H"), ]

print(dim(noDXvsDOX24H))
print(dim(noDXvsDOX42H))
print(dim(DOX24HvsDOX42H))

# Check genes
genes_noDoxvs24H <- noDXvsDOX24H[which(noDXvsDOX24H$Symbol %in% genes), c("Name", "Symbol", "Biotype", "DEG", "Direction", "log2FC", "padj")]
genes_noDoxvs42H <- noDXvsDOX42H[which(noDXvsDOX42H$Symbol %in% genes), c("Name", "Symbol", "Biotype", "DEG", "Direction", "log2FC", "padj")]
genes_DOX24HvsDOX42H <- DOX24HvsDOX42H[which(DOX24HvsDOX42H$Symbol %in% genes), c("Name", "Symbol", "Biotype", "DEG", "Direction", "log2FC", "padj")] # there are no significant gene in this comparison

# Save genes 
write.table(genes_noDoxvs24H, paste(dir_out, "/Check_genelist_AC19AC27_noDOXvs24H.txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(genes_noDoxvs42H, paste(dir_out, "/Check_genelist_AC19AC27_noDOXvs42H.txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)