library(readr)
library(TraMineR)
library(WeightedCluster)

df_school <- read_delim("Data/school.csv", delim = ";")
df_reprints <- read_delim("Data/reprints.csv", delim = ";")

# Clustering

## filtering data from 1919 to 2022
col_index <- which(colnames(df_reprints) == "1919")
school <- as.matrix(df_school[, col_index:214])
reprints <- as.matrix(df_reprints[, col_index:214])

## converting school data into sequence format
school.seq <- seqdef(school)

## calculating distances 
dschool <- seqdist(school.seq, method="OM", sm="INDELSLOG")

## measuring the quality of PAM clustering 
pamRange<-wcKMedRange(dschool, kvals = 2:10)
summary(pamRange, max.rank = 3)
plot(pamRange, stat = "all", norm = "zscore") # optimal number of clusters equals 4 (see fig. 5 in the paper)

## indentifying school and non-school clusters
school.clusters <- wcKMedoids(dschool, k = 4)
seqplot(school.seq, type='I', group = school.clusters$clustering, border = NA) # the biggest cluster "2026" consists of non-school novels

## clustering non-school novels
cluster_1_indices <- which(school.clusters$clustering == 2026)
reprints_cluster_1 <- reprints[cluster_1_indices,]

reprints.seq <- seqdef(reprints_cluster_1)
dreprints <- seqdist(reprints.seq, method="OM", sm="INDELSLOG")

pamRange2<-wcKMedRange(dreprints, kvals = 2:10)
summary(pamRange2, max.rank = 3)
plot(pamRange2, stat = "all", norm = "zscore") # optimal number of clusters equals 3 (see fig. 6 in the paper)
reprints.clusters <- wcKMedoids(dreprints, k = 3)

## merging school and and non-school clusters
final_clusters <- school.clusters$clustering
final_clusters[cluster_1_indices] <- reprints.clusters$clustering + 3
final_clusters <- factor(final_clusters, levels = c(147, 2008, 498, 560, 934, 996), labels = c("1. Occasional reprints", "2. Archive", "3. Soviet-Russian curricula", "4. Second canon", "5. Early Soviet curricula", "6. Soviet curricula"))

## visualising multidomain sequence plot
reprints.seq.full <- seqdef(reprints)
seqdom <- list("school domain"=school.seq, "reprints domain"=reprints.seq.full) 
require(colorspace) # opting for the palettes 
sch.col <- sequential_hcl(3, palette = "PurpOr")
rep.col <- sequential_hcl(3, palette = "Plasma")
color_pal=list(sch.col, rep.col)

dir.create("output", showWarnings = FALSE) 
png(filename = "output/clustering_results.png", width = 6000, height = 4500, units = "px", pointsize=65)
seqplotMD(seqdom, type="I", cpal.dom=color_pal, group = final_clusters, border = NA)
dev.off()


# Calculating integrative potential of each novel(1919-2022)

## integrative potential for reprints
seqindic.result_rep_all <- seqindic(reprints.seq.full, indic='integr', ipos.args=list(pos.states=c("1+","1 reprint")))
seqindic.result_rep_all <- cbind(df2[, 1:2], seqindic.result_rep_all)
colnames(seqindic.result_rep_all)[which(colnames(seqindic.result_rep_all) == "Integr")] <- "Integr_reprints"

## integrative potential for school
seqindic.result_school_all <- seqindic(school.seq, indic='integr', ipos.args=list(pos.states=c("curricula")))
seqindic.result_school_all <- cbind(df1[, 1:2], seqindic.result_school_all)
colnames(seqindic.result_school_all)[which(colnames(seqindic.result_school_all) == "Integr")] <- "Integr_school"

# plot the values of integrative potential in the reprints and school domains for each novel 
data <- left_join(seqindic.result_rep_all, seqindic.result_school_all, by = c("title", "author"))
ggplot(data, aes(x = Integr_reprints, y = Integr_school)) +
  geom_point(alpha = 0.5) +
  labs(x = "Integr reprints",
       y = "Integr school") +
  theme_linedraw() +
  theme(aspect.ratio = 1)


# Integrative potential per cluster (Soviet vs Russian post-Soviet period)

## Soviet period

col_1991 <- which(names(df_reprints) == "1991")
col_index <- which(names(df_reprints) == "1919")

reprints_Sov <- as.matrix(df_reprints[, col_index:col_1991])
reprints.seq_Sov <- seqdef(reprints_Sov)
seqindic.result_rep_Sov <- seqindic(reprints.seq_Sov, indic=c('integr'), ipos.args=list(pos.states=c("1+","1 reprint")))
seqindic.result_rep_Sov <- cbind(df_reprints[, 1:2], seqindic.result_rep_Sov)
colnames(seqindic.result_rep_Sov)[which(colnames(seqindic.result_rep_Sov) == "Integr")] <- "Integr_reprints_Sov"

school_Sov <- as.matrix(df_school[, col_index:col_1991])
school.seq_Sov <- seqdef(school_Sov)
seqindic.result_school_Sov <- seqindic(school.seq_Sov, indic=c('integr'), ipos.args=list(pos.states=c("curricula")))
seqindic.result_school_Sov <- cbind(df_school[, 1:2], seqindic.result_school_Sov)
colnames(seqindic.result_school_Sov)[which(colnames(seqindic.result_school_Sov) == "Integr")] <- "Integr_school_Sov"


## post-Soviet period
col_1992 <- which(names(df_reprints) == "1992")
col_2022 <- which(names(df_reprints) == "2022")

### reprint domain
reprints_postSov <- as.matrix(df_reprints[, col_1992:col_2022])
reprints.seq_postSov <- seqdef(reprints_postSov)
seqindic.result_rep_postSov <- seqindic(reprints.seq_postSov, indic=c('integr'), ipos.args=list(pos.states=c("1+","1 reprint")))
seqindic.result_rep_postSov <- cbind(df2[, 1:2], seqindic.result_rep_postSov)
colnames(seqindic.result_rep_postSov)[which(colnames(seqindic.result_rep_postSov) == "Integr")] <- "Integr_reprints_postSov"

### school domain
school_postSov <- as.matrix(df_school[, col_1992:col_2022])
school.seq_postSov <- seqdef(school_postSov)
seqindic.result_school_postSov <- seqindic(school.seq_postSov, indic=c('integr'), ipos.args=list(pos.states=c("curricula")))
seqindic.result_school_postSov <- cbind(df_school[, 1:2], seqindic.result_school_postSov)
colnames(seqindic.result_school_postSov)[which(colnames(seqindic.result_school_postSov) == "Integr")] <- "Integr_school_postSov"


## merging and visualising
d_integr_all <- data.frame(
  Integr_reprints_Sov = seqindic.result_rep_Sov$Integr_reprints_Sov,
  Integr_rep_postSov = seqindic.result_rep_postSov$Integr_reprints_postSov,
  Cluster = final_clusters
)
d_integr_all <- melt(d_integr_all, id.vars = "Cluster", variable.name = "Type", value.name = "Value")

### boxplot: Integrative potential of reprint sequences in each cluster over two periods: 1919-1991 and 1992-2022 (see fig. 4 in the paper)
palette <- brewer.pal(6, "Set2")
ggplot(d_integr_all, aes(x = Type, y = Value, fill = Cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = palette) +
  labs(x = "",
       y = "Value",
       title = "") +
  theme_linedraw()