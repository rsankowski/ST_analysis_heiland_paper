# load packages
library(tidyverse)
library(clusterProfiler)
library(readxl)
library(Seurat)
library(ggpubr)
library(SPATA2)

# ggplot2 plots can be easily combined with 'patchwork'
library(patchwork) 

# load objects

gbm_list <- map(file.path("objects",list.files("objects")),loadSpataObject)
gbm_list <- map(gbm_list,function(x) try(runDeAnalysis(x, across = "segmentation")))
map(gbm_list, function(x) try(sum(x@data[[1]]$scaled["LTBR",] > 0) ))
map(gbm_list, function(x) try(sum(x@data[[1]]$scaled["LTB",] > 0) ))


names(gbm_list) <- file.path("objects",list.files("objects"))
gbm_degene_list <- map(gbm_list, function(x) try(getDeaResultsDf(object = x, 
                                                              across = "segmentation", 
                                                              method_de = "wilcox")))
names(gbm_degene_list) <- file.path("objects",list.files("objects"))

gbm_degene_list <- gbm_degene_list[unlist(map(gbm_degene_list,is_tibble))] %>% 
  bind_rows(.id="sample")

human_gbm <- loadSpataObject("objects/spata_object_212_T.RDS")
human_gbm_2 <- loadSpataObject("objects/spata_object_275_T.RDS")
human_gbm_3 <- loadSpataObject("objects/spata_object_217_T.RDS")
human_gbm_4 <- loadSpataObject("objects/spata_object_218_T.RDS")

## add the LTB and LTBR status
gbm_list <- map(gbm_list, function(x) {
  tryCatch({
  .df <- getFeatureVariables(x, features = "segmentation", return = "data.frame")
  .df$LTB_status <- case_when(x@data[[1]]$counts["LTB",] > 0 ~ "LTB_positive",
                             T ~ "LTB_negative")
  .df$LTBR_status <- case_when(x@data[[1]]$counts["LTBR",] > 0 ~ "LTBR_positive",
                               T ~ "LTBR_negative")
  
  x <- addFeatures(x,.df,overwrite = T)
  
  x <- runDeAnalysis(x, across = "LTB_status")
  runDeAnalysis(x, across = "LTBR_status")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

## extract the diffgenes LTB spots
gbm_degene_list_ltb <- map(gbm_list, function(x) try(getDeaResultsDf(object = x, 
                                                 across = "LTB_status", 
                                                 method_de = "wilcox")))
names(gbm_degene_list_ltb) <- file.path("objects",list.files("objects"))

gbm_degene_list_ltb <- gbm_degene_list_ltb[unlist(map(gbm_degene_list_ltb,is_tibble))] %>% 
  bind_rows(.id = "sample") 
gbm_degene_list_ltb <- gbm_degene_list_ltb %>% 
  filter(p_val_adj<.05)

write.csv(gbm_degene_list_ltb,file.path("data","diffgenes_in_ltb_positive_spots.csv"))

##look at most common genes LTB spots
gbm_degene_list_ltb %>% 
  group_by(gene) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count)) %>% 
  write.csv(file.path("data","unique_diffgenes_in_ltb_positive_spots.csv"))

## GO term analysis of the genes in LTB spots
genes <- bitr(unique(gbm_degene_list_ltb$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
colnames(genes)[1] <- "gene"

## comparisob between datasets
gbm_degene_list_ltb_2 <- gbm_degene_list_ltb %>% 
  left_join(genes[!duplicated(genes$gene),]) %>% 
  na.omit()

ccomp <- compareCluster(ENTREZID ~ sample,
               data=gbm_degene_list_ltb_2, 
               fun = enrichGO,
               OrgDb = 'org.Hs.eg.db')
clusterProfiler::dotplot(ccomp)

enr_go <- ego <- enrichGO(gene          = gbm_degene_list_ltb_2$ENTREZID,
                          OrgDb         = 'org.Hs.eg.db',
                          minGSSize     = 5,
                          ont           = "BP",
                          pool          = TRUE,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
d <- as.data.frame(enr_go@result) %>% 
  filter(p.adjust<.05)

write.csv(d, file.path("data","go_term_analysis_LTBR_containing_dots.csv"))

d %>% 
  ggplot(aes(reorder(Description, -p.adjust), -log10(p.adjust), fill=-log10(p.adjust), size=Count)) +
  geom_point(color="black", stroke=.25, shape=21) +
  scale_fill_distiller(type="seq",palette = "Oranges", direction = 1) +
  coord_flip() +
  theme_light() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()) +
  expand_limits(y=c(1.2, 1.8))
ggsave(file.path("plots", "GO_terms_LTB_containing_dots.pdf"), height =4, width=11)


## LTBR
## extract the diffgenes
gbm_degene_list_ltbr <- map(gbm_list, function(x) try(getDeaResultsDf(object = x, 
                                                                      across = "LTBR_status", 
                                                                      method_de = "wilcox")))
names(gbm_degene_list_ltbr) <- file.path("objects",list.files("objects"))

gbm_degene_list_ltbr <- gbm_degene_list_ltbr[unlist(map(gbm_degene_list_ltbr,is_tibble))] %>% 
  bind_rows(.id = "sample") 
gbm_degene_list_ltbr <- gbm_degene_list_ltbr %>% 
  filter(p_val_adj<.05)

write.csv(gbm_degene_list_ltbr,file.path("data","diffgenes_in_ltbr_positive_spots.csv"))

##look at most common genes
gbm_degene_list_ltbr %>% 
  group_by(gene) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count)) %>% 
  write.csv(file.path("data","unique_diffgenes_in_ltbr_positive_spots.csv"))

## GO term analysis of the genes
genes <- bitr(unique(gbm_degene_list_ltbr$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
colnames(genes)[1] <- "gene"

## comparisob between datasets
gbm_degene_list_ltbr_2 <- gbm_degene_list_ltbr %>% 
  left_join(genes[!duplicated(genes$gene),]) %>% 
  na.omit()

ccomp <- compareCluster(ENTREZID ~ sample,
                        data=gbm_degene_list_ltbr_2, 
                        fun = enrichGO,
                        OrgDb = 'org.Hs.eg.db')
clusterProfiler::dotplot(ccomp)

enr_go <- ego <- enrichGO(gene          = gbm_degene_list_ltbr_2$ENTREZID,
                          OrgDb         = 'org.Hs.eg.db',
                          minGSSize     = 5,
                          ont           = "BP",
                          pool          = TRUE,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
d <- as.data.frame(enr_go@result) %>% 
  filter(p.adjust<.05)

write.csv(d, file.path("data","go_term_analysis_ltbrR_containing_dots.csv"))

d %>% 
  filter(grepl("(immune)", Description)) %>% 
  ggplot(aes(reorder(Description, -p.adjust), -log10(p.adjust), fill=-log10(p.adjust), size=Count)) +
  geom_point(color="black", stroke=.25, shape=21) +
  scale_fill_distiller(type="seq",palette = "Oranges", direction = 1) +
  coord_flip() +
  theme_light() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()) +
  expand_limits(y=c(1,2.6))
ggsave(file.path("plots", "filtered_GO_terms_ltbr_containing_dots.pdf"), height =3, width=8)
#human_gbm <- loadSpataObject("data/spata-obj-gbm275.RDS")

# store example genes of interest as character vectors
genes_a <- c("LTB", "LTBR")
genes_b <- c("CARTPT", "OLIG1", "GFAP", "SYNPR", "HOPX", "CCK")

plots <- plotSurfaceInteractive(object = human_gbm_2)

# compare gene expression on the surface
plotSurfaceComparison(object = human_gbm_4, 
                      color_by = genes_a,
                      smooth = TRUE, 
                      smooth_span = 0.2, 
                      pt_size = 3, 
                      pt_clrsp = "inferno")

## plot segmentation
plotSegmentation(human_gbm_2, pt_size = 2)

walk(names(gbm_list), function(x) {
  tryCatch({
  item <- gbm_list[[x]]
  pdf(file.path("plots",paste0(gsub(".*/","",x),"_segmeng_plot_LTB_LTBR.pdf")))
  print(plotSegmentation(item, pt_size = 3))
  dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  })

pdf(file.path("plots","segmeng_plots_LTB_LTBR.pdf"))
seg_plots
dev.off()

walk(names(gbm_list), function(x) {
  tryCatch({
    item <- gbm_list[[x]]
    pdf(file.path("plots",paste0(gsub(".*/","",x),"_surface_plot_LTB_LTBR.pdf")))
    print(plotSurfaceComparison(object = item, 
                                color_by = genes_a,
                                smooth = TRUE, 
                                smooth_span = 0.1, 
                                pt_size = 3, 
                                pt_clrsp = "inferno"))
    dev.off()
    
    jpeg(file.path("plots",paste0(gsub(".*/","",x),"_HE_image.jpg")))
    plot(gbm_list[[x]]@images[[1]])
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

walk(names(gbm_list), function(x) {
  tryCatch({
    item <- gbm_list[[x]]
    pdf(file.path("plots",paste0(gsub(".*/","",x),"_surface_plot_LTB_LTBR_image.pdf")))
    print(plotSurfaceComparison(object = item, 
                                color_by = genes_a,
                                pt_alpha = 0,
                                image =SPATA2::image(item),
                                smooth = TRUE, 
                                smooth_span = 0.1, 
                                pt_size = 3, 
                                pt_clrsp = "inferno"))
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

## load verhaak genesets
genesets <- map(list.files("/Users/romansankowski/Documents/single_cell_analysis/gl261_meta_analysis/data/MSIGDB_genesets", pattern = ".gmt$"), function(x) {read.gmt(file.path("/Users/romansankowski/Documents/single_cell_analysis/gl261_meta_analysis/data/MSIGDB_genesets",x))}) %>% 
  bind_rows() %>% 
  mutate(term=gsub("VERHAAK_GLIOBLASTOMA_","",term))


rownames(human_gbm@data[[1]]$counts[human_gbm@data[[1]]$counts>0,])

walk(names(gbm_list), function(x) {
  tryCatch({
    item <- gbm_list[[x]]
    input_list <- list(Classical = genesets$gene[genesets$term=="CLASSICAL"], 
                       Mesenchymal = genesets$gene[genesets$term=="MESENCHYMAL"], 
                       Neural = genesets$gene[genesets$term=="NEURAL"], 
                       Proneural = genesets$gene[genesets$term=="PRONEURAL"]
    )
    pdf(file.path("plots",paste0(gsub(".*/","",x),"_surface_plot_verhaak_signatures.pdf")))
    print(plotSurfaceAverage(object = item, 
                             color_by = input_list, 
                             smooth = TRUE, 
                             pt_size = 3, 
                             pt_clrsp = "inferno")
      )
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})
## plot neftel et al gene modules
gene_modules <- read_excel("/Users/romansankowski/Documents/single_cell_analysis/gl261_meta_analysis/data/Neftel_et_al/other_data/1-s2.0-S0092867419306877-mmc2.xlsx", skip = 4) %>% 
  pivot_longer(MES2:`G2/M`, names_to = "module", values_to = "gene") %>% 
  arrange(module) %>% 
  na.omit() #%>% 
#filter(!grepl("(G1/S|G2/M)", module)) %>% 
#mutate(module=gsub("[1-2]$", "", module)) #%>% 
#mutate(module=gsub("NPC|OPC", "NPC/OPC", module))

walk(names(gbm_list), function(x) {
  tryCatch({
    item <- gbm_list[[x]]
    input_list <- list("Mesenchymal-like 1" = gene_modules$gene[gene_modules$module=="MES1"], 
                       "Mesenchymal-like 2" = gene_modules$gene[gene_modules$module=="MES2"], 
                       "Astrocyte-like" = gene_modules$gene[gene_modules$module=="AC"], 
                       "NPC-like 1" = gene_modules$gene[gene_modules$module=="NPC1"], 
                       "NPC-like 2" = gene_modules$gene[gene_modules$module=="NPC2"], 
                       "OPC-like" = gene_modules$gene[gene_modules$module=="OPC"]
    )
    pdf(file.path("plots",paste0(gsub(".*/","",x),"_surface_plot_neftel_et_al_signatures.pdf")))
    print(plotSurfaceAverage(object = item, 
                             color_by = input_list, 
                             smooth = TRUE, 
                             pt_size = 3, 
                             pt_clrsp = "inferno")
    )
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

## plot gene modules
gene_modules <- read_excel("/Users/romansankowski/Documents/single_cell_analysis/gl261_meta_analysis/data/Neftel_et_al/other_data/1-s2.0-S0092867419306877-mmc2.xlsx", skip = 4) %>% 
  pivot_longer(MES2:`G2/M`, names_to = "module", values_to = "gene") %>% 
  arrange(module) %>% 
  na.omit() #%>% 
#filter(!grepl("(G1/S|G2/M)", module)) %>% 
#mutate(module=gsub("[1-2]$", "", module)) #%>% 
#mutate(module=gsub("NPC|OPC", "NPC/OPC", module))

## plot collapsed neftel et al gene modules
gene_modules <- read_excel("/Users/romansankowski/Documents/single_cell_analysis/gl261_meta_analysis/data/Neftel_et_al/other_data/1-s2.0-S0092867419306877-mmc2.xlsx", skip = 4) %>% 
pivot_longer(MES2:`G2/M`, names_to = "module", values_to = "gene") %>% 
  arrange(module) %>% 
  na.omit() %>% 
  #filter(!grepl("(G1/S|G2/M)", module)) %>% 
  mutate(module=gsub("[1-2]$", "", module)) #%>% 
  #mutate(module=gsub("NPC|OPC", "NPC/OPC", module))

walk(names(gbm_list), function(x) {
  tryCatch({
    item <- gbm_list[[x]]
    input_list <- list("Mesenchymal-like1" = gene_modules$gene[gene_modules$module=="MES"], 
                       "Astrocyte-like" = gene_modules$gene[gene_modules$module=="AC"], 
                       "NPC-like" = gene_modules$gene[gene_modules$module=="NPC"], 
                       "OPC-like" = gene_modules$gene[gene_modules$module=="OPC"]
    )
    pdf(file.path("plots",paste0(gsub(".*/","",x),"_surface_plot_neftel_et_al_collapsed_signatures.pdf")))
    print(plotSurfaceAverage(object = item, 
                             color_by = input_list, 
                             smooth = TRUE, 
                             pt_size = 3, 
                             pt_clrsp = "inferno")
    )
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

## score the spatial counts for glioblastoma signatures
## assign neftel tumor cell state 
gene_modules <- read_excel("/Users/romansankowski/Documents/single_cell_analysis/gl261_meta_analysis/data/Neftel_et_al/other_data/1-s2.0-S0092867419306877-mmc2.xlsx", skip = 4) %>% 
  pivot_longer(MES2:`G2/M`, names_to = "module", values_to = "gene") %>% 
  arrange(module) %>% 
  na.omit() %>% 
  #filter(!grepl("(G1/S|G2/M)", module)) %>% 
  mutate(module=gsub("[1-2]$", "", module)) %>% 
  mutate(module=gsub("NPC|OPC", "NPC/OPC", module))

## subtype and cell cycle scoring
lst <- list(
  #s.features = gene_modules$gene[gene_modules$module=="G1/S"], 
  g2m.features = gene_modules$gene[gene_modules$module=="G2/M"], 
  ac.features=gene_modules$gene[gene_modules$module=="AC"],
  mes.features=gene_modules$gene[gene_modules$module=="MES"],
  npc.opc.features=gene_modules$gene[gene_modules$module=="NPC/OPC"]
)

gbm_scores <- map(gbm_list, function(x) {
  tryCatch({
  .df <- x@data[[1]]$counts %>% 
      CreateSeuratObject() %>% 
      #SCTransform(variable.features.n = 10000) %>% 
      NormalizeData() %>% 
      ScaleData() 
  lst2 <- map(lst, function(x) x[x %in% rownames(df)]) 
  .df %>% 
    AddModuleScore(features = lst2, 
                   name = c("S_score", "G2M_Score", "AC_Score", "Mes_Score","NPC_OPC_Score"),
                   assay = "RNA",
                   ctrl = 50,
                   nbin = 10)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  })


## plot mean expression per segment
data.frame("segment"= human_gbm_2@fdata$`275_T`$segmentation,  "expression"=c(human_gbm_2@data$`275_T`$scaled["LTB",])) %>% 
  ggplot(aes(segment, expression, fill=segment, color=segment)) +
  geom_jitter() +
  stat_compare_means(comparisons = list(c("Cellular", "Infiltrative"), c("Cellular","none"), c("Cellular","Vascular_Hyper")))
