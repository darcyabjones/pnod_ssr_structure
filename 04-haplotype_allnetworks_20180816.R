library(pegas)
library(ape)
library(openxlsx)
library(data.table)

chaps <- read.csv("03-haplotypes/haplotypeshp_curated_20180112.csv")
chaps_core <- read.csv("03-haplotypes/haplotypescorehp_curated.csv")
chaps$gene <- as.character(chaps$gene)
chaps$seqid <- as.character(chaps$seqid)
chaps$new_hap <- as.character(chaps$new_hap)

View(chaps)
View(chaps_core)
chaps_core$gene <- as.character(chaps_core$gene)
chaps_core$seqid <- as.character(chaps_core$seqid)
chaps_core$new_hap <- as.character(chaps_core$new_hap)

tox1_core <- ape::read.FASTA("03-haplotypes/tox1coreseqs_24haps.fasta")
names(tox1_core) <- sapply(strsplit(names(tox1_core), " "), FUN = function(x) x[[1]])
toxa_core <- ape::read.FASTA("03-haplotypes/toxAcoreseqs3_25haps.fasta")
names(toxa_core) <- sapply(strsplit(names(toxa_core), " "), FUN = function(x) x[[1]])
tox3_core <- ape::read.FASTA("03-haplotypes/tox3coreseqs_16haps.fasta")
names(tox3_core) <- sapply(strsplit(names(tox3_core), " "), FUN = function(x) x[[1]])

tox1_haps_core <- with(chaps_core, chaps_core[gene == "Tox1",])
toxa_haps_core <- with(chaps_core, chaps_core[gene == "ToxA",])
tox3_haps_core <- with(chaps_core, chaps_core[gene == "Tox3",])

row.names(tox1_haps_core) <- tox1_haps_core$seqid
row.names(toxa_haps_core) <- toxa_haps_core$seqid
row.names(tox3_haps_core) <- tox3_haps_core$seqid
View(tox1_haps_core)
View(toxa_haps_core)
View(tox3_haps_core)

tox1_haps_core <- tox1_haps_core[names(tox1_core),]
toxa_haps_core <- toxa_haps_core[names(toxa_core),]
tox3_haps_core <- tox3_haps_core[names(tox3_core),]

# Here we make the SnTox3network

tox3_haplo <- pegas::haplotype(tox3_core)
tox3_haplo
attributes(tox3_haplo)
attr(tox3_haplo, "dimnames")[[1]] <- c(
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[1]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[2]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[3]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[4]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[5]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[6]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[7]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[8]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[9]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[10]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[11]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[12]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[13]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[14]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[15]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[16]], "new_hap"])
)



d3 <- read.dna("03-haplotypes/tox3coreseqs_16haps.fasta", format="fasta")
h3 <- haplotype(d3)
attributes(h3)
attr(h3, "dimnames")[[1]] <- c(
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[1]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[2]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[3]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[4]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[5]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[6]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[7]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[8]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[9]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[10]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[11]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[12]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[13]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[14]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[15]], "new_hap"]),
      unique(tox3_haps_core[attr(tox3_haplo, "index")[[16]], "new_hap"])
)
h3 <- sort(h3, what = "label")
net <- haploNet(h3)


ind.hap <- with(stack(setNames(attr(h3, "index"), rownames(h3))), table(hap=ind, pop=rownames(d3)[values]))
View(ind.hap)

svg(filename = "04-haplotype_networks/tox3_network_20180112.svg")
plot(net, size=attr(net, 'freq'), scale.ratio=0.9, cex = 0.8, pie=ind.hap)
legend(2, -2, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2, cex=0.5)
dev.off()

plot(net, size=attr(net, 'freq'), scale.ratio=0.9, cex = 0.8, pie=ind.hap)

legend(2, -2, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2, cex=0.5)

###################################################
# Here we make the SnTox1network
tox1_haplo <- pegas::haplotype(tox1_core)
tox1_haplo
attributes(tox1_haplo)
attr(tox1_haplo, "dimnames")[[1]] <- c(
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[1]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[2]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[3]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[4]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[5]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[6]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[7]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[8]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[9]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[10]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[11]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[12]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[13]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[14]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[15]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[16]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[17]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[18]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[19]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[20]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[21]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[22]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[23]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[24]], "new_hap"])
)

d1 <- read.dna("03-haplotypes/tox1coreseqs_24haps.fasta", format="fasta")
h1 <- haplotype(d1)
attributes(h1)
attr(h1, "dimnames")[[1]] <- c(
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[1]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[2]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[3]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[4]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[5]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[6]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[7]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[8]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[9]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[10]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[11]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[12]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[13]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[14]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[15]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[16]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[17]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[18]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[19]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[20]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[21]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[22]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[23]], "new_hap"]),
      unique(tox1_haps_core[attr(tox1_haplo, "index")[[24]], "new_hap"]))

h1 <- sort(h1, what = "label")
net1 <- haploNet(h1)


ind.hap1 <- with(stack(setNames(attr(h1, "index"), rownames(h1))), table(hap=ind, pop=rownames(d1)[values]))

svg(filename = "04-haplotype_networks/tox1_network1_20180112.svg")
plot(net1, size=attr(net1, 'freq'), scale.ratio=0.9, cex = 0.8, pie=ind.hap1)
legend(3, 7, colnames(ind.hap1), col=rainbow(ncol(ind.hap1)), pch=19, ncol=2, cex=0.5)
dev.off()

plot(net1, size=attr(net1, 'freq'), scale.ratio=0.9, cex = 0.8, pie=ind.hap1)

legend(3, 7, colnames(ind.hap1), col=rainbow(ncol(ind.hap1)), pch=19, ncol=2, cex=0.5)

##################################################
# Here we make the SnToxAnetwork
toxa_haplo <- pegas::haplotype(toxa_core)
toxa_haplo
attributes(toxa_haplo)
attr(toxa_haplo, "dimnames")[[1]] <- c(
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[1]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[2]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[3]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[4]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[5]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[6]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[7]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[8]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[9]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[10]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[11]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[12]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[13]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[14]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[15]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[16]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[17]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[18]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[19]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[20]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[21]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[22]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[23]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[24]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[25]], "new_hap"]))

da <- read.dna("03-haplotypes/toxAcoreseqs3_25haps.fasta", format="fasta")
ha <- haplotype(da)
attributes(ha)
attr(ha, "dimnames")[[1]] <- c(
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[1]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[2]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[3]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[4]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[5]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[6]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[7]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[8]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[9]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[10]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[11]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[12]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[13]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[14]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[15]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[16]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[17]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[18]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[19]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[20]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[21]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[22]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[23]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[24]], "new_hap"]),
      unique(toxa_haps_core[attr(toxa_haplo, "index")[[25]], "new_hap"]))

ha <- sort(ha, what = "label")
neta <- haploNet(ha)


ind.hapa <- with(stack(setNames(attr(ha, "index"), rownames(ha))), table(hap=ind, pop=rownames(da)[values]))

svg(filename = "04-haplotype_networks/toxa_network_20180112.svg")
plot(neta, size=attr(neta, 'freq'), scale.ratio=0.4, cex = 0.6, pie=ind.hapa)
legend(3, 13, colnames(ind.hapa), col=rainbow(ncol(ind.hapa)), pch=19, ncol=2, cex=0.5)
dev.off()

plot(neta, size=attr(neta, 'freq'), scale.ratio=0.6, cex = 0.6, pie=ind.hapa)

legend(3, 13, colnames(ind.hapa), col=rainbow(ncol(ind.hapa)), pch=19, ncol=2, cex=0.5)



