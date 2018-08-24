#find cor.test bw 3 temporal clusters and top-grown wheat lines
wheatnclusters <- read_csv("data/topgrownwheat_3clusters.csv")
View(wheatnclusters)


wheatnclusters$wheat3_avg <- rowSums(wheatnclusters[2:6])/5
wheatnclusters$wheat4_avg <- rowSums(wheatnclusters[7:11])/5
write.csv(wheatnclusters, "01-clustering/wheatnclusters.csv")


cor.test(wheatnclusters$cluster3, wheatnclusters$wheat3_avg)
cor.test(wheatnclusters$cluster4, wheatnclusters$wheat4_avg)
cor.test(wheatnclusters$cluster5, wheatnclusters$Mace_5)




