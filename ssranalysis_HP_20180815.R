library("readr")
library("tibble")
library("dplyr")
library("tidyr")
library("adegenet")
library("poppr")
library("openxlsx")
library("ggplot2")
library("colorspace") 
library("ape")
library("data.table")
library("plyr") 
library(svglite)
library(phangorn)
#############################
# Make phylogenetic tree of all isolates including satgo and outliers
# Start with the 30SSRHP20170803allna for the general pylogenetic tree 

ssr184_allna <- read.csv("data/30SSRHP20170803allnacorrected.csv", header = TRUE)
View(ssr184_allna)
dim(ssr184_allna)

ssr184_allna[ssr184_allna %in% c("-", ".")] <- NA
ssr184_allna$year <- factor(ssr184_allna$year)
ssr184_allna$year <- as.integer(ssr184_allna$year)
ssr184_allna_raw <- ssr184_allna
cols <- names(ssr184_allna)[!names(ssr184_allna) %in% c("location", "year", "species", "isolate", "fungicide", "loc", "note")]
ssr184_allna[, cols] <- round(ssr184_allna[, cols])
setDT(ssr184_allna)
ssr_markerdata <- fread("data/markers.csv")
invisible({
      ssr_markerdata[, `:=`(motif = strsplit(motif, "/"),
            fprimer = NULL,
            rprimer = NULL)]
})
ssr_markerdata
invisible(ssr184_allna[, SNOD8 := NULL])

col2keep_ssrallna <- c("location","year", "species", "isolate", "fungicide", "loc", "note")

tmp <- as.matrix(ssr184_allna[, !col2keep_ssrallna, with = FALSE])
View(tmp)
row.names(tmp) <- ssr184_allna[, isolate]


ssr184_allna_gd <- df2genind(
      tmp,
      ploidy = 1,
      type = "codom",
      pop = ssr184_allna[, loc],
      strata = ssr184_allna[, c("location","year", "species", "isolate", "fungicide", "loc", "note"), with = FALSE]
)
ssr184_allna_gd

ssr184_allna_gc <- as.genclone(ssr184_allna_gd)
ssr184_allna_gc

loc_order <- names(ssr184_allna_gc$loc.n.all)
setkey(ssr_markerdata, "name")

# fix_replen(ssr185_allna_gc, replen = ssr_markerdata[loc_order, motif_length])
ssr184_allna_gctree <- bruvo.boot(
      pop = ssr184_allna_gc,
      replen = ssr_markerdata[loc_order, motif_length],
      tree = "upgma",
      showtree = TRUE
)
outgroup_labels <- c("15FG039", "15FG46", "WAC11137", "SG1", "W1_1", "NFR")
ssrs_gc_tree <- root(ssr184_allna_gctree, outgroup = outgroup_labels, resolve.root = TRUE)

attributes(ssr184_allna_gctree)
ssr184_allna_gctree$tip.label
ssr184_allna_gctree$node.labels

write.tree(ssr184_allna_gctree, "01_phylogeny tree/ssr184_allna_gctree1.nwk")

write.tree(ssr184_allna_gctree, file="ssr184_allna_gcNewickTree.tre")
write.nexus(ssr184_allna_gctree, file="ssr184_allna_gcNexusTree.nex")
plot.phylo(ssr184_allna_gctree, type = "fan")

distance_matrix <- as.matrix(bruvo.dist(
      pop = ssr184_allna_gc,
      replen = ssr_markerdata[loc_order, motif_length]
))

write.xlsx(distance_matrix, "01_phylogeny tree/ssrs184allna_tree_distances.xlsx", row.names=TRUE)

distance_matrix[1:5, 1:5]
###############################
# Aus pop ssr analysis
ssrs <- readr::read_csv("data/20180301DJ-30SSR.csv")

id_cols <- c("location", "year", "species", "isolate", "fungicide", "loc", "lat", "lon", "region", "note")
ssr_cols <- ssrs %>% select(-one_of(id_cols)) %>% colnames()
ssrs[ssr_cols][is.na(ssrs[ssr_cols])] <- 0

ssrs_wa <- ssrs %>% filter(species == "stago" & !location %in% c("USA", "France"))

ssrs_wa
dim(ssrs_wa)
ssrs_wa_clonecor <- ssrs_wa %>% filter(!isolate %in% c("WAC13630", "WAC13528"))
dim(ssrs_wa_clonecor)


ssrs_loc_year <- ssrs_wa %>% group_by(location, year) %>% summarise(freq = n()) 

g <- ggplot(data = ssrs_loc_year, aes(location)) +
      geom_col(aes(y = freq, fill = year)) +
      scale_fill_distiller(type = "div", palette = "Spectral", values = c(1.0, 0.8, 0.0)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ylab("frequency")

ggsave(plot = g, filename = "01-clustering/nsamples_vs_location_bar_coloured_year_wa.pdf")
g

g <- ggplot(data = ssrs_loc_year, aes(as.character(year))) +
      geom_col(aes(y = freq, fill = location)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ylab("frequency") +
      xlab("year")

ggsave(plot = g, filename = "01-clustering/nsamples_vs_year_bar_coloured_location_wa.pdf")
g

############
# this section is for clonecorrect exploration

duplicates <- ssrs_wa[ssrs_wa %>% select(one_of(ssr_cols)) %>% duplicated() | ssrs_wa %>% select(one_of(ssr_cols)) %>% duplicated(fromLast = TRUE), ]

duplicates

duplicate_table <- as.data.frame(duplicates %>% select(one_of(ssr_cols)))
row.names(duplicate_table) <- duplicates$isolate


dist(duplicate_table)

tibble2genind <- function(df, pop = "location"){
      id_cols <- c("location", "year", "species", "isolate", "fungicide", "loc", "lat", "lon", "region", "note")
      tmp <- as.matrix(df %>% select(-one_of(id_cols)))
      row.names(tmp) <- pull(df, isolate)
      
      gd <- adegenet::df2genind(
            tmp,
            ploidy = 1,
            type = "codom",
            pop = pull(df, pop),
            strata = df %>% select(one_of(id_cols))
      )
      return(gd)
}

gd_wa <- tibble2genind(ssrs_wa)
gd_wa


adegenet::tab(gd_wa)[1:5, 1:9]

print(paste("number of loci = ", adegenet::nLoc(gd_wa)))
print(paste("number of individuals = ", adegenet::nInd(gd_wa)))

nAllele <- function(gen){
      return(length(colnames(adegenet::tab(gen))))
}
print(paste("number of alleles = ", nAllele(gd_wa)))

loci_stats <- poppr::locus_table(gd_wa)
write.csv(loci_stats, file = "01-clustering/locus_table.csv")
loci_stats


# Keyword arguments
# cutoff - numeric. A number from 0 to 1 defining the minimum number of differentiating samples.
# MAF	- numeric. A number from 0 to 1 defining the minimum minor allele frequency. This is passed as the thresh parameter of isPoly.
gd_wa <- poppr::informloci(gd_wa, cutoff = 2/adegenet::nInd(gd_wa), MAF = 0.01)

noclone <- poppr::clonecorrect(gd_wa, strata = NA)
noclone

clones <- setdiff(row.names(gd_wa@tab), row.names(noclone@tab))
clones

sort(as.matrix(dist(gd_wa))["WAC13630",])[1:3]

filter(ssrs, isolate %in% c("WAC13630", "WAC13528"))

gd_wa <- gd_wa[!row.names(gd_wa@tab) %in% c("16FG162", "16FG170")]
gd_wa

## PCA visualisation

gd_wa_pca <- ade4::dudi.pca(
      gd_wa,
      nf = nAllele(gd_wa) - 1,
      scannf = FALSE,
      scale = FALSE,
      center = TRUE
)

gd_wa_pca
dim(gd_wa_pca$tab)

pca_variance <- tibble::tibble(eigen_value = gd_wa_pca$eig)
pca_variance <- pca_variance %>% mutate(variance_explained = 100 * (eigen_value / sum(eigen_value, na.rm = TRUE)))

# to make PCA contribution figure
g <- ggplot(data=pca_variance, aes(x=seq_along(variance_explained), y=variance_explained)) + 
      geom_col() +
      xlab("Principal component") +
      ylab("% variance explained")

ggsave(plot = g, filename = "01-clustering/PCA_variance_explained.pdf")
g

gd_wa_pca_vectors <- gd_wa_pca$l1
gd_wa_pca_vectors$isolate <- row.names(gd_wa_pca_vectors)
gd_wa_pca_vectors <- tibble::as.tibble(gd_wa_pca_vectors) %>%
      rename_all(function(x){gsub("RS", "PC", x)})

metadata <- tibble::as.tibble(strata(gd_wa)) %>% 
      mutate_all(funs(as.character)) %>% 
      mutate_at(c("year", "fungicide", "lat", "lon"), funs(as.numeric))

gd_wa_pca_vectors <- dplyr::left_join(gd_wa_pca_vectors, metadata, by="isolate") %>%
      select(colnames(metadata), everything())

write.csv(gd_wa_pca_vectors, "01-clustering/wa_isolates_pca_values.csv")

# Just showing first 2 PCs and metadata.
gd_wa_pca_vectors[,1:12]

###############
# To make a scatter plot for PC1nPC2 by locations
g <- ggplot(data = gd_wa_pca_vectors, aes(x = PC1, y = PC2, colour = location, shape = location)) +
      ggtitle("PC1 vs PC2 coloured by isolate location") +
      geom_point() +
      scale_shape_manual(values = 1:length(unique(gd_wa_pca_vectors$location)) %% 25) +
      theme(legend.position="right") +
      guides(fill = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))

ggsave(plot = g, filename = "01-clustering/PC1_vs_PC2_coloured_location.pdf")
g

# To make a scatter plot for PC1nPC2 by years
g <- ggplot(data = gd_wa_pca_vectors, aes(x = PC1, y = PC2, color = year)) +
      ggtitle("PC1 vs PC2 coloured by isolate year") +
      geom_point() +
      scale_color_distiller(type = "div", palette = "Spectral", values = c(1.0, 0.8, 0.0)) +
      theme(legend.position="right")

ggsave(plot = g, filename = "01-clustering/PC1_vs_PC2_coloured_year.pdf")

g

g <- ggplot(data = gd_wa_pca_vectors, aes(x = PC1, y = PC2, colour = factor(year), shape = factor(year))) +
      ggtitle("PC1 vs PC2 coloured by year") +
      geom_point() +
      scale_shape_manual(values = 1:length(unique(gd_wa_pca_vectors$year)) %% 25) +
      theme(legend.position="right") +
      guides(fill = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
ggsave(plot = g, filename = "01-clustering/PC1_vs_PC2_colourednshaped_year.pdf")
g

##############################
## Clustering with snapclust.
# To make Figure2 for the ssr paper - 2 table charts by locations and years


clust_5 <- adegenet::snapclust(gd_wa, pop.ini = "ward", k = 5)

tmp <- as.data.frame(clust_5$group)
names(tmp)[1] <- "cluster_5"
View(tmp)
tmp$isolate <- row.names(tmp)

write.csv(tmp, "01-clustering/clustering_wa_5_groups.csv", row.names = FALSE)

compoplot(clust_5)

svg(filename = "01-clustering/compoplot_clust_5.svg", width = 6, height = 6)
compoplot(clust_5)
dev.off()


clust_5_location <- table(gd_wa@strata$location, clust_5$group)
colnames(clust_5_location) <- c("cluster1","cluster2", "cluster3", "cluster4", "cluster5")
table.value(clust_5_location, col.labels = 1:5)
table.value(clust_5_location, col.labels = c("cluster1","cluster2", "cluster3", "cluster4", "cluster5"))

clust_5_year <- table(gd_wa@strata$year, clust_5$group)
table.value(clust_5_year, col.labels = c("cluster1","cluster2", "cluster3", "cluster4", "cluster5"))


#################
# HP try to mahe a table with sorted years
clust_5_year_test <- table(gd_wa@strata$year, clust_5$group)
write.table(clust_5_year_test, "01-clustering/clust5_year.txt", row.names = TRUE)
clust_5_year_sorted <- read.table("01-clustering/clust5_year_sorted.txt", header = TRUE)
colnames(clust_5_year_sorted) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")

clust_5_year_sorted[,1]
rownames(clust_5_year_sorted)
table.value(clust_5_year_sorted)

#####################################

# Prepare for DAPC analysis for 5 cluster from MLGD

ssrs_wa_clonecor <- ssrs_wa %>% filter(!isolate %in% c("16FG162", "16FG170"))
dim(ssrs_wa_clonecor)
dim(tmp)
dim(ssrs_wa)

ssrs_wa5_clonecor <- merge(ssrs_wa_clonecor, tmp, by= "isolate")
dim(ssrs_wa5_clonecor)
unique(ssrs_wa5_clonecor$cluster_5)

colnames(ssrs_wa5_clonecor)
View(ssrs_wa5_clonecor)

write.csv(ssrs_wa5_clonecor, "01-clustering/ssrs_wa_clonecor_5grp.csv", row.names = FALSE)


tibble2genind_5grp <- function(df, pop = "cluster_5"){
      id_cols <- c("location", "year", "species", "isolate", "fungicide", "loc", "lat", "lon", "region", "note", "cluster_5")
      tmp <- as.matrix(df %>% select(-one_of(id_cols)))
      row.names(tmp) <- pull(df, isolate)
      
      gd <- adegenet::df2genind(
            tmp,
            ploidy = 1,
            type = "codom",
            pop = pull(df, pop),
            strata = df %>% select(one_of(id_cols))
      )
      return(gd)
}

gd_wa5 <- tibble2genind_5grp(ssrs_wa5_clonecor)
gd_wa5

adegenet::tab(gd_wa5)[1:5, 1:9]
head(adegenet::tab(gd_wa5))
dim(adegenet::tab(gd_wa5))
print(paste("number of loci = ", adegenet::nLoc(gd_wa5)))
print(paste("number of individuals = ", adegenet::nInd(gd_wa5)))

gd_wa5 <- poppr::informloci(gd_wa5, cutoff = 2/adegenet::nInd(gd_wa5), MAF = 0.01)
print(paste("number of loci = ", adegenet::nLoc(gd_wa5)))
print(paste("number of individuals = ", adegenet::nInd(gd_wa5)))

set.seed(0418)

xy1 <- xvalDapc(gd_wa5, gd_wa5@strata$cluster_5, training.set = 0.9, n.rep = 50, scale = FALSE, center = TRUE)

xy1

print(paste("Number of PCs Achieving Lowest MSE = ", xy1$`Number of PCs Achieving Lowest MSE`)) #10

set.seed(0418)
dapc_ssrs_wa5 <- dapc(gd_wa5, n.pca = 10, n.da = 4)
dapc_ssrs_wa5

summary(dapc_ssrs_wa5)
round(head(dapc_ssrs_wa5$posterior),3)

myPal <- colorRampPalette(c("blue","gold","red","green", "purple"))
# da1 & da2
svg(filename = "01-clustering/scatter-da1&2_Ausna5_20180418.svg", width = 6, height = 6)
scatter(dapc_ssrs_wa5, 1,2, scree.pca = FALSE,   legend = TRUE, col=transp(myPal(6)), scree.da=TRUE,cell=1.5, cex=2, bg="white",cstar=0)
dev.off()

# da1 & da3
svg(filename = "01-clustering/scatter-da1&3_Ausna5_20180418.svg", width = 6, height = 6)
scatter(dapc_ssrs_wa5, 1,3, scree.pca = FALSE,   legend = TRUE, col=transp(myPal(6)), scree.da=TRUE,cell=1.5, cex=2, bg="white",cstar=0)
dev.off()

# da1 & da4
svg(filename = "01-clustering/scatter-da1&4_Ausna5_20180418.svg", width = 6, height = 6)
scatter(dapc_ssrs_wa5, 1,4, scree.pca = FALSE,   legend = TRUE, col=transp(myPal(6)), scree.da=TRUE,cell=1.5, cex=2, bg="white",cstar=0)
dev.off()

# da2 & da3
svg(filename = "01-clustering/scatter-da2&3_Ausna5_20180418.svg", width = 6, height = 6)
scatter(dapc_ssrs_wa5, 2,3, scree.pca = FALSE,   legend = TRUE, col=transp(myPal(6)), scree.da=TRUE,cell=1.5, cex=2, bg="white",cstar=0)
dev.off()

# da2 & da4
svg(filename = "01-clustering/scatter-da2&4_Ausna5_20180418.svg", width = 6, height = 6)
scatter(dapc_ssrs_wa5, 2,4, scree.pca = FALSE,   legend = TRUE, col=transp(myPal(6)), scree.da=TRUE,cell=1.5, cex=2, bg="white",cstar=0)
dev.off()

# da3 & da4
svg(filename = "01-clustering/scatter-da3&4_Ausna5_20180309.svg", width = 6, height = 6)
scatter(dapc_ssrs_wa5, 3,4, scree.pca = FALSE,   legend = TRUE, col=transp(myPal(6)), scree.da=TRUE,cell=1.5, cex=2, bg="white",cstar=0)
dev.off()

######################
# Redo Table charts with DAPC grouping information/correction
# A table was prepared where infor of 5 group assignment based on DAPC analysis was added to other original infor

fiveDAPCgroup <- read_csv("data/hap&DAPCgroup_20180629.csv")
View(fiveDAPCgroup)

clust_5dapc_year <- table(fiveDAPCgroup$Year, fiveDAPCgroup$DAPCgroup)

clust_5dapc_year <- as.tibble(clust_5dapc_year)
colnames(clust_5dapc_year) <- c("year", "group", "count")

labs <- list(xlab("Group"), ylab("Year"))

g <- ggplot(data = clust_5dapc_year, aes(x=group, y=year, size=count)) +
      geom_point(shape=15, na.rm = TRUE) +
      scale_size_continuous(name = "Count", limits = c(1, 30), breaks = c(1, 5, 10, 15, 20, 25)) +
      labs

ggsave("01-clustering/cluster_summary_year.svg", width = 2.8, height = 4.)
print(g)

clust_5dapc_location <- table(fiveDAPCgroup$Location, fiveDAPCgroup$DAPCgroup)

clust_5dapc_location <- as.tibble(clust_5dapc_location)
colnames(clust_5dapc_location) <- c("location", "group", "count")

labs <- list(xlab("Group"), ylab("Location"))

g <- ggplot(data = clust_5dapc_location, aes(x=group, y=location, size=count)) +
      geom_point(shape=15, na.rm = TRUE) +
      scale_size_continuous(name = "Count", limits = c(1, 20), breaks = c(1, 5, 10, 15, 20)) +
      labs

ggsave("01-clustering/cluster_summary_location.svg", width = 3.5, height = 5)

print(g)

######################
#getting infor on loading of significant alleles for LD1toLD4

loadingplot(dapc_ssrs_wa5$var.contr, axis=1, thres=0.03)
loadingplot(dapc_ssrs_wa5$var.contr, axis=2, thres=0.03)
loadingplot(dapc_ssrs_wa5$var.contr, axis=3, thres=0.03)
loadingplot(dapc_ssrs_wa5$var.contr, axis=4, thres=0.03)

loading_LD1to4 <- dapc_ssrs_wa5$var.contr
class(loading_LD1to4)

colnames(loading_LD1to4)

loading_LD1 <- loading_LD1to4[loading_LD1to4[,"LD1"] >= 0.03,]
loading_LD2 <- loading_LD1to4[loading_LD1to4[,"LD2"] >= 0.03,]
loading_LD3 <- loading_LD1to4[loading_LD1to4[,"LD3"] >= 0.03,]
loading_LD4 <- loading_LD1to4[loading_LD1to4[,"LD4"] >= 0.03,]
loadingsignificant_Ld1n2 <- rbind(loading_LD1, loading_LD2)
loadingsignificant_Ld1to3 <- rbind(loadingsignificant_Ld1n2, loading_LD3)
loadingsignificant_Ld1to4 <- rbind(loadingsignificant_Ld1to3, loading_LD4)
write.csv(loadingsignificant_Ld1to4, file = "01-clustering/loadingsignificant_Ld1to4.csv")


freq_ssr7 <- tab(genind2genpop(gd_wa5[loc=c("SSR7")]), freq = TRUE)
freq_ssr7[,colSums(freq_ssr7) != 0]
write.csv(freq_ssr7, file = "01-clustering/freq_ssr7.csv")

freq_ssr8 <- tab(genind2genpop(gd_wa5[loc=c("SSR8")]), freq = TRUE)
freq_ssr8[,colSums(freq_ssr8) != 0]
write.csv(freq_ssr8, file = "01-clustering/freq_ssr8.csv")

freq_ssr9 <- tab(genind2genpop(gd_wa5[loc=c("SSR9")]), freq = TRUE)
freq_ssr9[,colSums(freq_ssr9) != 0]
write.csv(freq_ssr9, file = "01-clustering/freq_ssr9.csv")

freq_ssr10 <- tab(genind2genpop(gd_wa5[loc=c("SSR10")]), freq = TRUE)
freq_ssr8[,colSums(freq_ssr8) != 0]
write.csv(freq_ssr10, file = "01-clustering/freq_ssr10.csv")

freq_ssr12 <- tab(genind2genpop(gd_wa5[loc=c("SSR12")]), freq = TRUE)
freq_ssr12[,colSums(freq_ssr12) != 0]
write.csv(freq_ssr12, file = "01-clustering/freq_ssr12.csv")

freq_ssr22 <- tab(genind2genpop(gd_wa5[loc=c("SSR22")]), freq = TRUE)
freq_ssr22[,colSums(freq_ssr22) != 0]
write.csv(freq_ssr22, file = "01-clustering/freq_ssr22.csv")

freq_ssr23 <- tab(genind2genpop(gd_wa5[loc=c("SSR23")]), freq = TRUE)
write.csv(freq_ssr23, file = "01-clustering/freq_ssr23.csv")

freq_ssr27 <- tab(genind2genpop(gd_wa5[loc=c("SSR27")]), freq = TRUE)
write.csv(freq_ssr27, file = "01-clustering/freq_ssr27.csv")

freq_ssr15 <- tab(genind2genpop(gd_wa5[loc=c("SSR15")]), freq = TRUE)
write.csv(freq_ssr15, file = "01-clustering/freq_ssr15.csv")

freq_ssr14 <- tab(genind2genpop(gd_wa5[loc=c("SSR14")]), freq = TRUE)
write.csv(freq_ssr14, file = "01-clustering/freq_ssr14.csv")
freq_ssr1 <- tab(genind2genpop(gd_wa5[loc=c("SSR1")]), freq = TRUE)
write.csv(freq_ssr1, file = "01-clustering/freq_ssr1.csv")
freq_ssr4 <- tab(genind2genpop(gd_wa5[loc=c("SSR4")]), freq = TRUE)
write.csv(freq_ssr4, file = "01-clustering/freq_ssr4.csv")
freq_ssr5 <- tab(genind2genpop(gd_wa5[loc=c("SSR5")]), freq = TRUE)
write.csv(freq_ssr5, file = "01-clustering/freq_ssr5.csv")
freq_ssr18 <- tab(genind2genpop(gd_wa5[loc=c("SSR18")]), freq = TRUE)
write.csv(freq_ssr18, file = "01-clustering/freq_ssr18.csv")
freq_ssr20 <- tab(genind2genpop(gd_wa5r[loc=c("SSR20")]), freq = TRUE)
write.csv(freq_ssr20, file = "01-clustering/freq_ssr20.csv")
freq_ssr22 <- tab(genind2genpop(gd_wa5[loc=c("SSR22")]), freq = TRUE)
write.csv(freq_ssr22, file = "01-clustering/freq_ssr22.csv")
freq_ssr28 <- tab(genind2genpop(gd_wa5[loc=c("SSR28")]), freq = TRUE)
write.csv(freq_ssr28, file = "01-clustering/freq_ssr28.csv")
freq_ssr19 <- tab(genind2genpop(gd_wa5[loc=c("SSR19")]), freq = TRUE)
write.csv(freq_ssr19, file = "01-clustering/freq_ssr19.csv")
freq_ssr25 <- tab(genind2genpop(gd_wa5[loc=c("SSR25")]), freq = TRUE)
write.csv(freq_ssr25, file = "01-clustering/freq_ssr25.csv")
freq_ssr27 <- tab(genind2genpop(gd_wa5[loc=c("SSR27")]), freq = TRUE)
write.csv(freq_ssr27, file = "01-clustering/freq_ssr27.csv")

###########################
# MAking 2 structure plots (from clust_5 and dapc) and information on dapc membership assignments
# And identify which isolates get different assignments with dapc compared to the snapclust() method
compoplot(clust_5)

clust_5_proba <- clust_5$proba
write.csv(clust_5_proba, file = "01-clustering/compoplot_proba_clust_5.csv")

svg(filename = "01-clustering/compoplot_clust_5.svg", width = 6, height = 6)
compoplot(clust_5)
dev.off()


dapc_ssrs_wa5_asignind <- assignplot(dapc_ssrs_wa5)
assignplot(dapc_ssrs_wa5)
assignplot(dapc_ssrs_wa5, subset = 1:50)
assignplot(dapc_ssrs_wa5, subset = 51:100)
assignplot(dapc_ssrs_wa5, subset = 101:153)


compoplot(dapc_ssrs_wa5, txt.leg=paste("Group", 1:5), ncol=1, col=c("blue", "green", "red", "purple", "yellow"),cleg=0.7)

svg(filename = "01-clustering/compoplot_ausdapc_5.svg", width = 6, height = 6)
compoplot(ssrs_wa5_clonecor_gccor_dapc, txt.leg=paste("Cluster", 1:5), ncol=1, col=c("blue", "green", "red", "purple", "yellow"),cleg=0.7)
dev.off()


####################
# Calculate amova table and Ia values for the whole pop and 5 groups

gc_gd_wa5 <- as.genclone(gd_wa5)
poppr(gc_gd_wa5, ~cluster_5)

table(strata(ssrs_wa5_clonecor_gccor, ~cluster_5))
gc_gd_wa5_amova <- poppr.amova(gc_gd_wa5, ~cluster_5)
gc_gd_wa5_amova

set.seed(0331)
gc_gd_wa5_amova_signiftest <- randtest(gc_gd_wa5_amova, nrepet = 999)
plot(gc_gd_wa5_amova_signiftest)
gc_gd_wa5_amova_signiftest


ia(gc_gd_wa5, sample = 999)

gc_gd_wa5_grp1 <- popsub(gc_gd_wa5, "1")
ia(gc_gd_wa5_grp1, sample = 999)

gc_gd_wa5_grp2 <- popsub(gc_gd_wa5, "2")
ia(gc_gd_wa5_grp2, sample = 999)

gc_gd_wa5_grp3 <- popsub(gc_gd_wa5, "3")
ia(gc_gd_wa5_grp3, sample = 999)

gc_gd_wa5_grp4 <- popsub(gc_gd_wa5, "4")
ia(gc_gd_wa5_grp4, sample = 999)

gc_gd_wa5_grp5 <- popsub(gc_gd_wa5, "5")
ia(gc_gd_wa5_grp5, sample = 999)


##########################################################
#To make weather plot

library(ggplot2) # must be version > 2.2.0

# import data
weather <- read_csv("data/weatherdata.csv")
# show first observations
head(weather)
weather$Avearge_Rainfall2 <- weather$Avearge_Rainfall/5
View(weather)


weatherplot <- ggplot(weather,aes(year)) +
      geom_line(aes(y= Min_Temperature, colour = "Min_Temperature"), linetype = "dashed") +
      geom_line(aes(y= Max_Temperature, colour = "Max_Temperature"), linetype = "dashed") +
      geom_line(aes(y= Avearge_Rainfall/5, colour = "Avearge_Rainfall"),linetype="dashed")+
      geom_smooth(method = "lm", se=FALSE, color="orange", aes(y = Min_Temperature)) +
      geom_smooth(method = "lm", se=FALSE, color="black", aes(y = Max_Temperature)) +
      geom_smooth(method = "lm", se=FALSE, color="red", aes(y = Avearge_Rainfall2))


# now adding the secondary axis, following the example in the help file ?scale_y_continuous
# and, very important, reverting the above transformation
weatherplot <- weatherplot + scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Avearge_Rainfall"))   

# modifying colours and theme options
weatherplot <- weatherplot + scale_colour_manual(values = c("blue", "purple", "dark green"))
weatherplot <- weatherplot + labs(y = "Temperature [°C]", x = "Year", colour = "Parameter")

weatherplot


yearNmaxtemp <- lm(weather$Max_Temperature ~ weather$year)

summary(yearNmaxtemp)
# y = -24.41 + 0.02x
summary(lm(weather$Min_Temperature ~ weather$year))
# y = 48.98 - 0.02x
summary(lm(weather$Avearge_Rainfall ~ weather$year))
# y = 594.59 - 0.27x


##########################################################
# pairwise genetic distance among 5 Aus groups 

gd_wa5_gp <- genind2genpop(gd_wa5, pop = ~cluster_5)

gd_wa5_gp_tree <- aboot(gd_wa5_gp, sample = 1000, quiet = TRUE, showtree = TRUE)

gd_wa5_gp_tree

write.tree(gd_wa5_gp_tree, "01-clustering/ssrs_wa5_tree.nwk")

gd_wa5_gp_distancematrix <- as.matrix(nei.dist(gd_wa5_gp))

write.xlsx(gd_wa5_gp_distancematrix, "01-clustering/ssrs_wa5_tree_distances_nei.xlsx", row.names=TRUE)

gd_wa5_gp_distancematrix_edwards <- as.matrix(edwards.dist(gd_wa5_gp))
write.xlsx(gd_wa5_gp_distancematrix_edwards, "01-clustering/ssrs_wa5_tree_distances_edwards.xlsx", row.names=TRUE)

gd_wa5_gp_distancematrix_rogers <- as.matrix(rogers.dist(gd_wa5_gp))
write.xlsx(ssr_wa5_distance_matrix_rogers, "01-clustering/ssrs_wa5_tree_distances_rogers.xlsx", row.names=TRUE)

gd_wa5_gp_distancematrix_reynolds <- as.matrix(reynolds.dist(gd_wa5_gp))
write.xlsx(ssr_wa5_distance_matrix_reynolds, "01-clustering/ssrs_wa5_tree_distances_raynolds.xlsx", row.names=TRUE)

gd_wa5_gp_distancematrix_provesti <- as.matrix(provesti.dist(gd_wa5_gp))
write.xlsx(ssr_wa5_distance_matrix_provesti, "01-clustering/ssrs_wa5_tree_distances_provesti.xlsx", row.names=TRUE)

# NOTE: "provesti" distance was chosen following the review paper Grunwald et al 2017 "Best practices for population analysis"

###########################################



