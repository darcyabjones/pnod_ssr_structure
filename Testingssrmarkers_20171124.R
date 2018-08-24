library("openxlsx")
library("ggplot2")
library("colorspace") # Used in clustergrams
library("mclust")
library("ape")
library("data.table")
library("plyr") # Used in clustergrams
library("adegenet")
library("svglite")
library(magrittr)
library(dplyr)
library(treemap)
library(mmod)
library(poppr)

comellelein11isolates <- read.csv("commonallelein11isolates_20171124.csv")
View(comellelein11isolates)
head(comellelein11isolates)
dim(comellelein11isolates)
WAC8410 <- as.numeric(comellelein11isolates[1,6:35])
Sn15 <- as.numeric(comellelein11isolates[2,6:35])
WAC13071 <- as.numeric(comellelein11isolates[3,6:35])
WAC13073 <- as.numeric(comellelein11isolates[4,6:35])
WAC13532 <- as.numeric(comellelein11isolates[5,6:35])
WAC13524 <- as.numeric(comellelein11isolates[6,6:35])
WAC13617 <- as.numeric(comellelein11isolates[7,6:35])
FG114 <- as.numeric(comellelein11isolates[8,6:35])
FG119 <- as.numeric(comellelein11isolates[9,6:35])
FG49 <- as.numeric(comellelein11isolates[10,6:35])
RAC2182 <- as.numeric(comellelein11isolates[11,6:35])
isolates_wps <- comellelein11isolates


isolates_wps$year <- factor(as.integer(isolates_wps$year))
isolates_wps[isolates_wps %in% c("-", ".")] <- NA
isolates_wps[is.na(isolates_wps)] <- 0
View(isolates_wps)

cols <- names(isolates_wps)[!names(isolates_wps) %in% c("isolate", "loc","Pop", "Subpop", "year")]
isolates_wps[, cols] <- lapply(isolates_wps[cols], as.character)
setDT(isolates_wps)
setkey(isolates_wps, isolate)

isolates_wps_mat <- as.matrix(isolates_wps[, cols, with = FALSE])
row.names(isolates_wps_mat) <- isolates_wps[, isolate]

isolates_wps_gd <- df2genind(
      isolates_wps_mat,
      ploidy = 1,
      type = "codom",
      ind.names = isolates_wps[, isolate],
      pop = isolates_wps[, Pop],
      strata = isolates_wps[,c("Pop", "Subpop")])

isolates_wps_gdcor <- informloci(isolates_wps_gd)
isolates_wps_gdcor <- clonecorrect(isolates_wps_gdcor, strata = NA)

nLoc(isolates_wps_gdcor) #28
nInd(isolates_wps_gdcor) #11

isolates_wps_gc <- as.genclone(isolates_wps_gdcor)
isolates_wps_gc


library(ape)
library(phangorn)


ssr_markerdata <- fread("data/markers.csv")
invisible({
      ssr_markerdata[, `:=`(motif = strsplit(motif, "/"),
            fprimer = NULL,
            rprimer = NULL)]
})
ssr_markerdata

loc_order <- names(isolates_wps_gc$loc.n.all)
setkey(ssr_markerdata, "name")


isolates_wps_gctree <- bruvo.boot(
      pop = isolates_wps_gc,
      replen = ssr_markerdata[loc_order, motif_length],
      tree = "upgma",
      showtree = TRUE
)

##########################################
library("openxlsx")
library("ggplot2")
library("colorspace") # Used in clustergrams
library("mclust")
library("ape")
library("data.table")
library("plyr") # Used in clustergrams
library("adegenet")
library("svglite")
library(magrittr)
library(dplyr)
library(treemap)
library(mmod)
library(poppr)
###########################
ssrAusnatime3_raw <- read.csv("data/Ausnatime3_20171124.csv", header = TRUE)
View(ssrAusnatime3_raw)
unique(ssrAusnatime3_raw$loc)
ssrAusnatime3 <- ssrAusnatime3_raw
dim(ssrAusnatime3) #[1] 155  35

ssrAusnatime3_ssr1$year <- factor(as.integer(ssrAusnatime3_ssr1$year))
ssrAusnatime3_ssr1[ssrAusnatime3_ssr1 %in% c("-", ".")] <- NA
ssrAusnatime3[is.na(ssrAusnatime3)] <- 0
ssrAusnatime3_raw <- ssrAusnatime3
View(ssrAusnatime3)

###########
# Ausnatime3 without SSR1 marker
ssrAusnatime3_ssr1 <- ssrAusnatime3[,c(1:5, 7:35)]
View(ssrAusnatime3_ssr1)

cols <- names(ssrAusnatime3_ssr1)[!names(ssrAusnatime3_ssr1) %in% c("isolate", "loc","Pop", "Subpop", "year")]
ssrAusnatime3_ssr1[, cols] <- lapply(ssrAusnatime3_ssr1[cols], as.character)
setDT(ssrAusnatime3_ssr1)
setkey(ssrAusnatime3_ssr1, isolate)

ssrAusnatime3ssr1_mat <- as.matrix(ssrAusnatime3_ssr1[, cols, with = FALSE])
row.names(ssrAusnatime3ssr1_mat) <- ssrAusnatime3_ssr1[, isolate]

ssrAusnatime3ssr1_gd <- df2genind(
      ssrAusnatime3ssr1_mat,
      ploidy = 1,
      type = "codom",
      ind.names = ssrAusnatime3_ssr1[, isolate],
      pop = ssrAusnatime3_ssr1[, Pop],
      strata = ssrAusnatime3_ssr1[,c("Pop", "Subpop")])

ssrAusnatime3ssr1_gdcor <- informloci(ssrAusnatime3ssr1_gd)
ssrAusnatime3ssr1_gdcor <- clonecorrect(ssrAusnatime3ssr1_gdcor, strata = NA)
nLoc(ssrAusnatime3ssr1_gdcor) #27
nInd(ssrAusnatime3ssr1_gdcor) #152

strata(ssrAusnatime3ssr1_gdcor) <- strata(ssrAusnatime3ssr1_gdcor)[, c(1:2)]
popsummary_Ausnatime3_ssr1 <- poppr(ssrAusnatime3ssr1_gdcor)
write.csv(popsummary_Ausnatime3_ssr1, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr1_20171125.csv")

##############
# Ausnatime3 without SSR2 marker
ssrAusnatime3_ssr3 <- ssrAusnatime3[,c(1:6, 8:35)]
View(ssrAusnatime3_ssr3)

cols <- names(ssrAusnatime3_ssr3)[!names(ssrAusnatime3_ssr3) %in% c("isolate", "loc","Pop", "Subpop", "year")]
ssrAusnatime3_ssr3[, cols] <- lapply(ssrAusnatime3_ssr3[cols], as.character)
setDT(ssrAusnatime3_ssr3)
setkey(ssrAusnatime3_ssr3, isolate)

ssrAusnatime3ssr3_mat <- as.matrix(ssrAusnatime3_ssr3[, cols, with = FALSE])
row.names(ssrAusnatime3ssr3_mat) <- ssrAusnatime3_ssr3[, isolate]

ssrAusnatime3ssr3_gd <- df2genind(
      ssrAusnatime3ssr3_mat,
      ploidy = 1,
      type = "codom",
      ind.names = ssrAusnatime3_ssr3[, isolate],
      pop = ssrAusnatime3_ssr3[, Pop],
      strata = ssrAusnatime3_ssr3[,c("Pop", "Subpop")])

ssrAusnatime3ssr3_gdcor <- informloci(ssrAusnatime3ssr3_gd)
ssrAusnatime3ssr3_gdcor <- clonecorrect(ssrAusnatime3ssr3_gdcor, strata = NA)
nLoc(ssrAusnatime3ssr3_gdcor) #27
nInd(ssrAusnatime3ssr3_gdcor) #152

strata(ssrAusnatime3ssr3_gdcor) <- strata(ssrAusnatime3ssr3_gdcor)[, c(1:2)]
popsummary_Ausnatime3_ssr3 <- poppr(ssrAusnatime3ssr3_gdcor)
write.csv(popsummary_Ausnatime3_ssr3, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr3_20171125.csv")

##############
# Ausnatime3 without SSR4 marker
ssrAusnatime3_ssr4 <- ssrAusnatime3[,c(1:7, 9:35)]
View(ssrAusnatime3_ssr4)

write.csv(popsummary_Ausnatime3_ssr4, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr4_20171125.csv")

#############
# Ausnatime3 without SSR5 marker
ssrAusnatime3_ssr5 <- ssrAusnatime3[,c(1:8, 10:35)]
View(ssrAusnatime3_ssr5)

write.csv(popsummary_Ausnatime3_ssr5, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr5_20171125.csv")

##############
# Ausnatime3 without SSR6 marker
ssrAusnatime3_ssr6 <- ssrAusnatime3[,c(1:9, 11:35)]
View(ssrAusnatime3_ssr6)


write.csv(popsummary_Ausnatime3_ssr6, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr6_20171125.csv")

###############
# Ausnatime3 without SSR7 marker
ssrAusnatime3_ssr7 <- ssrAusnatime3[,c(1:10, 12:35)]
View(ssrAusnatime3_ssr7)

write.csv(popsummary_Ausnatime3_ssr7, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr7_20171125.csv")

#############
# Ausnatime3 without SSR8 marker
ssrAusnatime3_ssr8 <- ssrAusnatime3[,c(1:11, 13:35)]
View(ssrAusnatime3_ssr8)

write.csv(popsummary_Ausnatime3_ssr8, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr8_20171125.csv")

#############
# Ausnatime3 without SSR9 marker
ssrAusnatime3_ssr9 <- ssrAusnatime3[,c(1:12, 14:35)]
View(ssrAusnatime3_ssr9)

write.csv(popsummary_Ausnatime3_ssr9, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr9_20171125.csv")

###########
# Ausnatime3 without SSR10 marker
ssrAusnatime3_ssr10 <- ssrAusnatime3[,c(1:13, 15:35)]
View(ssrAusnatime3_ssr10)

write.csv(popsummary_Ausnatime3_ssr10, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr10_20171125.csv")

#########
# Ausnatime3 without SSR12 marker
ssrAusnatime3_ssr12 <- ssrAusnatime3[,c(1:14, 16:35)]
View(ssrAusnatime3_ssr12)

write.csv(popsummary_Ausnatime3_ssr12, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr12_20171125.csv")

##########
# Ausnatime3 without SSR18 marker
ssrAusnatime3_ssr18 <- ssrAusnatime3[,c(1:15, 17:35)]
View(ssrAusnatime3_ssr18)

write.csv(popsummary_Ausnatime3_ssr18, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr18_20171125.csv")

############
# Ausnatime3 without SSR20 marker
ssrAusnatime3_ssr20 <- ssrAusnatime3[,c(1:16, 18:35)]
View(ssrAusnatime3_ssr20)

write.csv(popsummary_Ausnatime3_ssr20, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr20_20171125.csv")

############
# Ausnatime3 without SSR21 marker
ssrAusnatime3_ssr21 <- ssrAusnatime3[,c(1:17, 19:35)]
View(ssrAusnatime3_ssr21)

write.csv(popsummary_Ausnatime3_ssr21, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr21_20171125.csv")

############
# Ausnatime3 without SSR22 marker
ssrAusnatime3_ssr22 <- ssrAusnatime3[,c(1:18, 20:35)]
View(ssrAusnatime3_ssr22)

write.csv(popsummary_Ausnatime3_ssr22, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr22_20171125.csv")

############
# Ausnatime3 without SSR23 marker
ssrAusnatime3_ssr23 <- ssrAusnatime3[,c(1:19, 21:35)]
View(ssrAusnatime3_ssr23)

write.csv(popsummary_Ausnatime3_ssr23, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr23_20171125.csv")

############
# Ausnatime3 without SSR26 marker
ssrAusnatime3_ssr26 <- ssrAusnatime3[,c(1:20, 22:35)]
View(ssrAusnatime3_ssr26)

write.csv(popsummary_Ausnatime3_ssr26, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr26_20171125.csv")

############
# Ausnatime3 without SSR27 marker
ssrAusnatime3_ssr27 <- ssrAusnatime3[,c(1:21, 23:35)]
View(ssrAusnatime3_ssr27)

write.csv(popsummary_Ausnatime3_ssr27, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr27_20171125.csv")

############
# Ausnatime3 without x marker
ssrAusnatime3_ssr28 <- ssrAusnatime3[,1:34]
View(ssrAusnatime3_ssr28)

write.csv(popsummary_Ausnatime3_ssr28, "18-Ausnatime3_withouteachssr_20171125/popsummary_Ausnatime3ssr24_20171125.csv")
# This was done until all 30 markers were taken, only chnage the 1st line of col chosen and the name of the output
#############################

library(magrittr)


allpopsumary_files <- list.files(".", pattern="popsummary_Ausnatime3*") %>% lapply(read.csv)

popsummary_all <- do.call(rbind, allpopsumary_files)
dim(popsummary_all)
# [1] 120  14
write.csv(popsummary_all, "./popsummary_all.csv")
View(popsummary_all)
popsummary_all2 <- popsummary_all[,c(2:3,7:9, 11:13)]
View(popsummary_all2)
write.csv(popsummary_all2, "./popsummary_all2.csv")

View(ssrAusnatime3)
markernames <- names(ssrAusnatime3)[6:35]

popsumary_1972to1996 <- popsummary_all2[popsummary_all2$Pop == "1972to1996",]
head(popsumary_1972to1996)
dim(popsumary_1972to1996) #[1] 30  8
View(popsumary_1972to1996)
rownames(popsumary_1972to1996) <- markernames
write.csv(popsumary_1972to1996, "./popsumary_1972to1996.csv")

popsumary_2001to2013 <- popsummary_all2[popsummary_all2$Pop == "2001to2013",]
head(popsumary_2001to2013)
dim(popsumary_2001to2013) #[1] 30  8
rownames(popsumary_2001to2013) <- markernames
write.csv(popsumary_2001to2013, "./popsumary_2001to2013.csv")

popsumary_2014to2016 <- popsummary_all2[popsummary_all2$Pop == "2014to2016",]
head(popsumary_2014to2016)
dim(popsumary_2014to2016) #[1] 30  8
rownames(popsumary_2014to2016) <- markernames
write.csv(popsumary_2014to2016, "./popsumary_2014to2016.csv")

popsumary_Total <- popsummary_all2[popsummary_all2$Pop == "Total",]
head(popsumary_Total)
dim(popsumary_Total) #[1] 30  8
rownames(popsumary_Total) <- markernames
write.csv(popsumary_Total, "./popsumary_Total.csv")

hist(popsumary_Total$Hexp)
totalmorethan715 <- popsumary_Total[popsumary_Total$Hexp > 0.715,]
dim(totalmorethan715) #[1] 6 8
write.csv(totalmorethan715, "./totalmorethan715.csv")

totallessthan715 <- popsumary_Total[popsumary_Total$Hexp < 0.715,]
write.csv(totalmorethan715, "./totalmorethan715.csv")
dim(totallessthan715) # [1] 24  8

shapiro.test(popsumary_Total$Hexp)
qqnorm(popsumary_Total$Hexp);qqline(popsumary_Total$Hexp, col = 2)

hist(popsumary_Total$H)
hist(popsumary_Total$G)
hist(popsumary_Total$lambda)
hist(popsumary_Total$Ia)
hist(popsumary_Total$rbarD)
