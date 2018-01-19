##### Aanderud Lab Phyloseq code ######
######### Metagenome Analysis #########
## Last updated: 11/13/17 ##
## Compiled by EFJ from http://joey711.github.io/phyloseq-demo/phyloseq-demo.html and google


############### Download and load packages  ###############
# source('http://bioconductor.org/biocLite.R')     ##For downloading phyloseq, only needs to be done once (regular install not available for phyloseq)
# biocLite('phyloseq')                             ##For downloading phyloseq, only needs to be done once (regular install not available for phyloseq)
library("phyloseq")
library("ggplot2")
library("ggthemes")
library("vegan")
library(MASS)
library("cluster")
require("RColorBrewer")
library("grid")
library("cooccur")
library(igraph)
library(Hmisc)
library(minerva)
library(parallel)

setwd("~/Box Sync/BYU/metagenomics/mothurproducts/102017")  ## Set to current working directory
theme_set(theme_bw())                              ## I also like theme_calc, theme_classic, theme_gray


################ Load and curate data ####################
community = import_mothur(mothur_list_file = NULL, mothur_group_file = NULL,
                  mothur_tree_file = NULL, cutoff = NULL,
                  "./conn.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared",
                  "./conn.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy")
                                                                 ### Read in mothur files, can also take tre file
print(community)                                                 ### Check summary of newly created phyloseq object
colnames(tax_table(community))<-c( "Domain", "Phylum",  "Class", 
                                "Order", "Family",  "Genus" )     ### Fix table to have correct taxa labels
rank_names(community)                                             ### check that it worked
get_taxa_unique(community, "Class")                               ### List observed taxa for a designated level
#community=subset_taxa(community, Class!= "Chloroplast")           ### Remove chloroplast taxa (leaves cyano and cyano_chloroplast_unclassified)
design <- read.delim("./conn_design102017.txt", row.names=1)      ### read in design .txt file
community<-merge_phyloseq(community, sample_data(design))         ### merge sample data with phyloseq object
sample_variables(community)                                       ### List of variable names from design file
sample_names(community)                                           ### List of sample names from design file
community_rare<- rarefy_even_depth(community, 
        sample.size = min(sample_sums(community)),rngseed = 711)          ##Rarefy data by smallest sample size
communityrel=transform_sample_counts(community, function(x) 100*x/sum(x))  ### Relativize data to relative abundance
community_rare_rel<- transform_sample_counts(community_rare, function(x) 100*x/sum(x))



############### Calculate and Plot Richness ##################
### Additional examples at http://joey711.github.io/phyloseq/plot_richness-examples
richness=estimate_richness(community_rare, measures=c("Observed","Shannon"))       ### Calculate richness values for observed and shannon
write.csv(richness, file="richness.csv")                                  ##Export to csv
write.csv(design, file="design.csv")                                  ##Export to csv
richness=read.csv("richness.csv", row.names = "X")                                    ### Reimport csv 
design2=read.csv("design.csv", row.names = "X")                             ### I have to do this because when I try to merge them originally, I lose the top sample and I don't know why
richnesstest=merge(richness,design2, by="row.names" )                         ### Merge richness with sample metadata in design file
aov1 = aov (Shannon~Watershed*Loc*Month_name, data=richnesstest)              ### ANOVA on Shannon data
summary(aov1)
tukey1 =  TukeyHSD(aov1, conf.level=0.95)
tukey1


aov2 = aov (Observed~Watershed*Loc*Month_name, data=richnesstest)             ### ANOVA on Observed data
summary(aov2)
tukey2 =  TukeyHSD(aov2, conf.level=0.95)
tukey2

## Graph using phyloseq to show individual points
# plot_richness(community_rare, x= "Watershed", color="Month_name", shape="Location", measures=c("Observed","Shannon"))              ### Plot individual sample richness by variable x
# plot_richness(community_rare, x= "Watershed", color="Loc", shape="Watershed", measures=c("Observed","Shannon"))
# plot_richness(community_rare, x= "Loc", color="Watershed", shape="Effort", measures=c("Observed","Shannon"))

## Graph using ggplot to make bargraphs
locationcolor=c("green3","blue2","red3")                                                      ### Select colors equal to number of unique values for "fill" in code below
qplot(Loc, Shannon, fill=Watershed, data=richnesstest, geom="boxplot", 
      ylab="Shannon Diversity",xlab="Location") +theme_bw() + scale_fill_manual(values=locationcolor)                          ### Plot richness data in box and whisker 
qplot(Loc, Observed, fill=Watershed, data=richnesstest, geom="boxplot", 
      ylab="Observed Richness",xlab="Location") +theme_bw() + scale_fill_manual(values=locationcolor)                   ### Plot richness data in box and whisker 



############ Network distance ##############
### Additional examples at http://joey711.github.io/phyloseq/plot_network-examples
set.seed(711L)                                                  ### Set random number generator for distance plot
plot_net(community_rare, maxdist = .6, color = "Loc", shape="Watershed")                    ### Plot network--- sometimes this doesn't work and I don't know why, use below instead
plot_network(make_network(community_rare, dist.fun="bray", max.dist = 0.65), community_rare, color = "Loc", shape="Watershed")



################ Heatmap ###################
### Additional examples at http://joey711.github.io/phyloseq/plot_heatmap-examples
wh0 = genefilter_sample(community_rare_rel, filterfun_sample(function(x) x > 0.005), A=0.5*nsamples(community_rare_rel))
community1 = prune_taxa(wh0, community_rare_rel)
communitypt <- prune_taxa(names(sort(taxa_sums(community_rare),TRUE)[1:100]), community_rare)    ### subset data to top XXX taxa [1:XXX]
### Make heatmap of subset data using phyloseq
# plot_heatmap(community_rare_rel, "Effort", taxa.label = "Phylum")                               
plot_heatmap(communitypt, "PCoA", "bray", "Effort", taxa.label = "Order")

### Make heat map using basic stats graphing instead of phyloseq - includes dendrogram
glom <- tax_glom(community_rare_rel, taxrank="Family")         ### Make object summarized by taxonomic rank
write.csv((otu_table(glom)), file='otusper-Family.csv')        ### Write csv of otu table
write.csv((tax_table(glom)), file='taxper-Family.csv')         ### Write csv of tax table
     ## Replace taxonomy names from taxper into otusper file and trim out low abundance OTUS in Excel, save as heatmap-class ##
     ## You can also average columns on treatments 
heatmapdata <- read.csv("heatmap-Family-avg.csv")                  ### Reimport merged data
rnames <- heatmapdata[,1]                                     ### assign labels in column 1 to "rnames"
mat_data <- data.matrix(heatmapdata[,2:ncol(heatmapdata)])    ### transform columns into a matrix
rownames(mat_data) <- rnames                                  ### assign row names 
my_palette <- colorRampPalette(c("white","red4"))(n = 200)    ### Make color ramp for heatmap graph
heatmap(mat_data, col=my_palette)                             ### Graph heatmap



############## Distance Statistics (ADONIS) ##############
set.seed(1)                                               ### set random number generator for repeatability
community_bray <- phyloseq::distance(community_rare_rel, method = "bray")    ### Calculate bray curtis distance matrixc
bad=c("13_2_d" , "FD_11_d" , "HS_2_d" , "HS_5_d")
community_rare_relnoFD <- subset_samples(community_rare_rel, !(DNACode %in%  bad))
community_braynoFD <- phyloseq::distance(community_rare_relnoFD, method = "bray")    ### Calculate bray curtis distance matrixc

communitydf <- data.frame(sample_data(community_braynoFD))                   ### make a data frame from the sample_data
adonis(community_bray ~ Watershed*Loc*Month_name, data = communitydf)   ### Adonis test
mod <- with(communitydf, betadisper(community_bray, Watershed))
mod

############### PCOA ###################
### http://joey711.github.io/phyloseq/plot_ordination-examples
community.ord = ordinate(community_rare_rel, "PCoA", "bray", weighted=TRUE)                     ###Ordinate phyloseq object
ord=plot_ordination(community_rare_rel, community.ord, color="Location", shape="Watershed")     ### Make plot
#ord+geom_polygon(aes(fill=Location, alpha=0.5)) + geom_point(size=3) + ggtitle("Communities-bray")  ### Plot polygons
ord + geom_vline(aes(xintercept = 0), lty=3, lwd=.5 )+ geom_hline(aes(yintercept = 0), lty=3, lwd=.5 )+ geom_point(size=2.25)+   theme_few() +
  scale_color_manual(values=c("forestgreen","springgreen3", "red", "steelblue1","blue"))        ### Plot manual color ramp

csin=hclust(community_bray, method="single")
caver=hclust(community_bray, method="aver")
ccom=hclust(community_bray, method="complete")
 plot(ccom, hang = -1)
plot(caver, hang = -1)
plot(csin, hang = -1)
rect.hclust(csin, 4)
plot(ccom, hang = -1)
 rect.hclust(ccom, 5)
plot(caver, hang = -1)
 rect.hclust(caver, 4)
 cl <- cutree(ccom, 4)

 ord <- cmdscale(community_bray)
 ordiplot(ord) 
 ordihull(ord, cl)
 
 
 fig<- ordiplot(ord, dis = "si")
 text(fig, "sites", cex=0.6)
ordihull(ord, cutree(ccom, 2))

ordiplot(ord, dis = "si")
 ordicluster(ord, ccom)
 

############ RDA or CCA ###############

## varetrim is a design file with environmental data in it -- must all be numeric and there cannot be empty cells, header in column 1 should be "Samples"
#varetrim=read.csv("varetrimnoFD.csv",  row.names = "Samples", stringsAsFactors = FALSE)      ## first version
#varetrim=read.csv("varetrimnoFD_2.csv",  row.names = "Samples", stringsAsFactors = FALSE)    ## to run isoMDS
varetrim=read.csv("varetrimnoFD_sigvect.csv", row.names = "Samples", stringsAsFactors = FALSE)   ## Only significant factors from isoMDS
vare.mds= isoMDS(community_braynoFD)                                 ### calculate NMDS
ef <- envfit(vare.mds, varetrim, permu = 999)                    ### calculate fit of NMDS
ef                                                               ### View statistics on envfit

## Edit otusper-class or heatmap-class file to have desired sample names and header name in column 1 save as vare-xxx##
vareotu= read.csv("heatmap-FamilynoFD.csv",  row.names = "Family")
vareotu<-as.data.frame(t(vareotu), stringsAsFactors = FALSE)     ## Transpose otu table to match format

##Calculate RDA
myrda <- rda(vareotu, varetrim)                                 ### Replace rda with cca to run CCA
myrda
anova(myrda)
## Plot RDA ##
plot(myrda, type="n")                                    ### open new plot
points(myrda, col=rdacolors, pch=rdashapes, cex=1.25)              ### Add points for samples (can also use text())
points(myrda, dis="cn",  cex=1.25)                      ### Add vectors for environmental data
text(myrda, "species", col="blue", cex=0.5)              ### Add names of OTUs (can also use points())

#Calculate CCA
mycca <- cca(vareotu, varetrim)                                 ### Replace rda with cca to run CCA
mycca
anova(mycca)
## Plot cca ##
plot(mycca, type="n", xlim=c(-3,2))                                    ### open new plot
points(mycca, col=rdacolors, pch=rdashapes, cex=1.25)              ### Add points for samples (can also use text())
points(mycca, dis="cn", cex=1.5)                      ### Add vectors for environmental data
text(mycca, "species", col="blue", cex=0.5)              ### Add names of OTUs (can also use points())
envfit(mycca,design,permutations = 999)

#color and shape exclude FD11
rdacolors= c("springgreen3","springgreen3","springgreen3","red","red","steelblue1","steelblue1","steelblue1","blue","blue","blue","forestgreen","forestgreen","blue","blue","red","forestgreen","forestgreen","forestgreen","steelblue1","steelblue1","steelblue1","blue","blue","blue","steelblue1","steelblue1","red","red","red","forestgreen","forestgreen","forestgreen","springgreen3","springgreen3","springgreen3","red","red")
rdashapes= c( 15 , 15 , 15 , 17 , 17 , 15 , 15 , 15 , 17 , 17 , 17 , 16 , 16 , 15 , 15 , 17 , 15 , 15 , 15 , 17 , 17 , 17 , 16 , 16 , 16 , 16 , 16 , 15 , 15 , 15 , 17 , 17 , 17 , 16 , 16 , 16 , 16 , 16 )

varemetal=read.csv("varemetal.csv", row.names = "Samples", stringsAsFactors = FALSE)   ## Only significant factors from isoMDS
vare.mds= isoMDS(community_braynoFD)                                 ### calculate NMDS
ef <- envfit(vare.mds, varemetal, permu = 999)                    ### calculate fit of NMDS
ef                                                               ### View statistics on envfit

## Edit otusper-class or heatmap-class file to have desired sample names and header name in column 1 save as vare-xxx##
vareotu= read.csv("otus_metals.csv",  row.names = "Family")
vareotu<-as.data.frame(t(vareotu), stringsAsFactors = FALSE)     ## Transpose otu table to match format

##Calculate RDA
myrda <- rda(vareotu, varemetal)                                 ### Replace rda with cca to run CCA
myrda
anova(myrda)
## Plot RDA ##
plot(myrda, type="n")                                    ### open new plot
points(myrda, col=rdacolors, pch=rdashapes, cex=1.25)              ### Add points for samples (can also use text())
points(myrda, dis="cn",  cex=1.25)                      ### Add vectors for environmental data
text(myrda, "species", col="blue", cex=0.5)              ### Add names of OTUs (can also use points())


############# Abundance plots ###############
#### Additional examples at http://joey711.github.io/phyloseq/plot_bar-examples

community.chl=subset_taxa(community, Phylum== "Chloroplast")  ### Make subset of data with only select taxa to check chloroplast and cyano abundance
community.cyo=subset_taxa(community, Class=="Cyanobacteria")
cyopt<-subset_taxa(community, Class=="Cyanobacteria")
bar_rel=transform_sample_counts(cyopt, function(x) 100*x/sum(x))     ### Relativize pruned data
mdf = psmelt(bar_rel)                                                     ### New dataframe of OTU table
mdf = psmelt(subset_samples(bar_rel,DNACode !="13_2_d"))

#### Relativized abundances ####
# plot_bar(community_rare, "Phylum")                                       ### Make abundance graph with relativized data

communitypt <- prune_taxa(names(sort(taxa_sums(community_rare),TRUE)[1:50]), community_rare)  ## Prune to top 25 taxa
bar_rel=transform_sample_counts(communitypt, function(x) 100*x/sum(x))     ### Relativize pruned data
mdf = psmelt(bar_rel)                                                     ### New dataframe of OTU table
mdf = psmelt(subset_samples(bar_rel,DNACode !="13_2_d"))                ### New dataframe of OTU table (Remove extra RB sample)

### Plot relativized abundance using ggplot instead of phyloseq 
p = ggplot(mdf, aes(x=Location, y=Abundance, fill=Genus, color=Genus))
p + geom_bar( stat="identity", position="stack") + facet_grid(Season~Watershed) + 
  scale_fill_manual(values=tol42rainbow)+ scale_color_manual(values=tol42rainbow)   

###Use whatever rainbow color brewer is large enough to cover the number of unique taxa
tol42rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",
                "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
 tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
 tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tolsiterainbow=c("#E8601C","#E8601C","#E8601C","#1965B0","#1965B0","#F6C141","#F6C141","#F6C141","#7BAFDE","#7BAFDE","#7BAFDE","#117744","#117744","#F7EE55","#F7EE55","#1965B0","#DC050C","#DC050C","#DC050C","#5289C7","#5289C7","#5289C7","#CAE0AB","#CAE0AB","#CAE0AB","#90C987","#90C987","#F1932D","#F1932D","#F1932D","#114477","#114477","#114477","#4EB265","#4EB265","#4EB265","#44AA77","#44AA77") 
# Function for plotting colors side-by-side
pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
pal(tolsiterainbow)
pal(tol18rainbow)





############### NMDS with ellipses ###################
otu_mat<-t(mat_data)                                        ### transpose mat_data
NMDS.log<-log1p(otu_mat)                                    ### take log
set.seed(42)                                                ### set random number generator for repeatability
sol <- metaMDS(otu_mat)                                     ### calculate MDS
NMDS <- data.frame(MDS1 = sol$points[,1], 
                   MDS2 = sol$points[,2],group=design$Loc)  ### Turn MDS into a data frame, where "group" is unique identifier for replicates
plot(sol$points, col = design$Loc)                          ### Start new plot
ord<-ordiellipse(sol, design$Loc, display = "sites", 
                 kind = "se", conf = 0.95, label = T)       ### Calculate ellipse on "groups"
df_ell <- data.frame()                                      ### make new df
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
          vegan:::veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}                                                           ### Merge and format NMDS into new df
ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=2)   ### Plot NMDS with geom_path for ellipses




############## Cooccurance ###############
# # Extract abundance matrix from the phyloseq object
# OTU1 = as(otu_table(community), "matrix")
# # transpose if necessary
# if(taxa_are_rows(physeq1)){OTU1 <- t(OTU1)}
# # Coerce to data.frame
# OTUdf = as.data.frame(OTU1)
# cooccur_mat<- data.matrix(glom)

OTUdf<- otu_table(communitypt)
dataPA <-(OTUdf > 0)*1 
dataPA<- as.data.frame(dataPA, row.names="X")
conn.cooccur <- cooccur(dataPA, type="spp_site")  ### Presence/absence abundance data frame or matrix
class(conn.cooccur)

# Return a Summary of Significant Results
summary(conn.cooccur)

# Return a Table of Results
conn.cooccur.tab <- prob.table(conn.cooccur)

# Plot a Visual Represnetation of Species Co-Occurrences
plot(conn.cooccur)

# If We are Interested in a 'Focal Specie", We can focus on that with this code:
pair(mod=conn.cooccur.analysis, spp="")
pair(mod=conn.cooccur.analysis, spp="")



############### Trees ###################
### http://joey711.github.io/phyloseq/plot_tree-examples



# ########### Gap statistic ##########
# ### I don't actually know what any of this means, but the code works and makes a graph :)
# pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}
# x = phyloseq:::scores.pcoa(community.ord, display="sites")
# gskmn = clusGap(x[, 1:2], FUN=pam1, K.max = 6, B = 50)
# gskmn
# 
# plot_clusgap = function(clusgap, title="Gap Statistic calculation results"){
#   gstab = data.frame(clusgap$Tab, k=1:nrow(clusgap$Tab))
#   p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
#   p = p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
#   p = p + ggtitle(title)
#   return(p)
# }
# gap_statistic_ordination = function(ord, FUNcluster, type="sites", K.max=6, axes=c(1:2), B=500, verbose=interactive(), ...){
#   require("cluster")
#   #   If "pam1" was chosen, use this internally defined call to pam
#   if(FUNcluster == "pam1"){
#     FUNcluster = function(x,k) list(cluster = pam(x, k, cluster.only=TRUE))     
#   }
#   # Use the scores function to get the ordination coordinates
#   x = phyloseq:::scores.pcoa(ord, display=type)
#   #   If axes not explicitly defined (NULL), then use all of them
#   if(is.null(axes)){axes = 1:ncol(x)}
#   #   Finally, perform, and return, the gap statistic calculation using cluster::clusGap  
#   clusGap(x[, axes], FUN=FUNcluster, K.max=K.max, B=B, verbose=verbose, ...)
# }
# gs = gap_statistic_ordination(community.ord, "pam1", B=50, verbose=FALSE)
# plot_clusgap(gs)






######### GEPHI Networks ###########
## Subset communities by watershed
communityRB<-subset_samples(community, Watershed=="Red Butte")
communityLR<-subset_samples(community, Watershed=="Logan")
communityPR<-subset_samples(community, Watershed=="Provo")

## Subset watershed by location
above<-subset_samples(communityRB, Loc=="Above")
above_rare<- rarefy_even_depth(above, 
                               sample.size = min(sample_sums(above)),rngseed = 711)          ##Rarefy data by smallest sample size

urban<-subset_samples(communityRB, Loc=="Urban")
urban_rare<- rarefy_even_depth(urban, 
                               sample.size = min(sample_sums(urban)),rngseed = 711)          ##Rarefy data by smallest sample size
dam<-subset_samples(communityRB, Loc=="Dam")
dam_rare<- rarefy_even_depth(dam, 
                             sample.size = min(sample_sums(dam)),rngseed = 711)          ##Rarefy data by smallest sample size
write.table(otu_table(above_rare), file="otu_RBabove", sep="\t")
write.table(otu_table(urban_rare), file="otu_RBurban", sep="\t")
write.table(otu_table(dam_rare), file="otu_RBdam", sep="\t")

above<-subset_samples(communityPR, Loc=="Above")
above_rare<- rarefy_even_depth(above, 
                               sample.size = min(sample_sums(above)),rngseed = 711)          ##Rarefy data by smallest sample size

urban<-subset_samples(communityPR, Loc=="Urban")
urban_rare<- rarefy_even_depth(urban, 
                               sample.size = min(sample_sums(urban)),rngseed = 711)          ##Rarefy data by smallest sample size
dam<-subset_samples(communityPR, Loc=="Dam")
dam_rare<- rarefy_even_depth(dam, 
                             sample.size = min(sample_sums(dam)),rngseed = 711)          ##Rarefy data by smallest sample size

write.table(otu_table(above_rare), file="otu_PRabove", sep="\t")
write.table(otu_table(urban_rare), file="otu_PRurban", sep="\t")
write.table(otu_table(dam_rare), file="otu_PRdam", sep="\t")

above<-subset_samples(communityLR, Loc=="Above")
above_rare<- rarefy_even_depth(above, 
                               sample.size = min(sample_sums(above)),rngseed = 711)          ##Rarefy data by smallest sample size

urban<-subset_samples(communityLR, Loc=="Urban")
urban_rare<- rarefy_even_depth(urban, 
                               sample.size = min(sample_sums(urban)),rngseed = 711)          ##Rarefy data by smallest sample size
dam<-subset_samples(communityLR, Loc=="Dam")
dam_rare<- rarefy_even_depth(dam, 
                             sample.size = min(sample_sums(dam)),rngseed = 711)          ##Rarefy data by smallest sample size

write.table(otu_table(above_rare), file="otu_LRabove", sep="\t")
write.table(otu_table(urban_rare), file="otu_LRurban", sep="\t")
write.table(otu_table(dam_rare), file="otu_LRdam", sep="\t")



#Load data
data<-read.table("otu_RBdam",header=T,row.names=1, sep="\t")   ## no
data<-read.table("otu_PRdam",header=T,row.names=1, sep="\t")   ## yes
data<-read.table("otu_LRdam",header=T,row.names=1, sep="\t")   ## no
data<-read.table("otu_RBurban",header=T,row.names=1, sep="\t") ## yes
data<-read.table("otu_LRurban",header=T,row.names=1, sep="\t") ## no
data<-read.table("otu_PRurban",header=T,row.names=1, sep="\t") ## yes
data<-read.table("otu_LRabove",header=T,row.names=1, sep="\t") ## yes
data<-read.table("otu_PRabove",header=T,row.names=1, sep="\t") ## yes
data<-read.table("otu_RBabove",header=T,row.names=1, sep="\t") ## yes
data=t(data) #transposes data to meet specifications

#head(data) #returns just first few rows of data
str(data) #describes the types of data in your rmatrix
#data


#filter (taxa designation for the number of taxa present in samples 
#20%=.20 50%=.50, 1=100% of the sites)
min.occu<-.75
data.filt<-data[, apply(data==0, 2, sum) < nrow(data)-min.occu*nrow(data)]

#transform from 0 to 1
range.0to1<- function(x){(x-min(x))/(max(x)-min(x))}
data.range<-apply(data.filt,2,range.0to1)
#head(data.range)

data.mine <- mine(data.range)

#reformat results
data.mine$MIC[upper.tri(data.mine$MIC, diag=T)] <- NA

data.results<-as.data.frame(cbind(which(!is.na(data.mine$MIC), arr.ind=T, useNames=F), na.omit(as.vector(data.mine$MIC))))
data.results[,1] <- rownames(data.mine$MIC)[data.results[,1]]
data.results[,2] <- rownames(data.mine$MIC)[data.results[,2]]
colnames(data.results) <- c("x", "y", "MIC")

#filter results by MIC threshold
mic <- 0.7
data.results.filt <-subset(data.results, data.results[,3] > mic)

#function to transform to a graph object
graph.transform.weights <- function (X) {
  require(igraph)
  data.tmp <- matrix(0, nrow=nrow(X), ncol=2)
  dimnames(data.tmp)[[2]] <- c("i", "j")
  data.tmp[,1] <- X[,1]
  data.tmp[,2] <- X[,2]
  data <- data.frame(data.tmp)
  graph.data <- graph.data.frame(data, directed=F)
  E(graph.data)$weight <- abs(as.numeric(X[,3]))
  E(graph.data)$rho<- as.numeric(X[,3])
  summary(graph.data)
  cat("Average degree:",ecount(graph.data)/vcount(graph.data)*2)
  return(graph.data)
}

data.graph<-graph.transform.weights(data.results.filt)

#write graph for gephi
#write.graph(data.graph, "data_network_RBdam.graphml", format="graphml")
#write.graph(data.graph, "data_network_PRdam.graphml", format="graphml")
#write.graph(data.graph, "data_network_LRdam.graphml", format="graphml")
#write.graph(data.graph, "data_network_RBurban.graphml", format="graphml")
#write.graph(data.graph, "data_network_LRurban.graphml", format="graphml")
#write.graph(data.graph, "data_network_PRurban.graphml", format="graphml")
#write.graph(data.graph, "data_network_LRabove.graphml", format="graphml")
#write.graph(data.graph, "data_network_PRabove.graphml", format="graphml")
#write.graph(data.graph, "data_network_RBabove.graphml", format="graphml")
write.graph

#################
##I stoppped here
#################


############
#Pearsons###
############
#Load data
data<-read.table("20140224_antarctica_one_at_a_time.txt",header=TRUE,sep="\t")
head(data_raw) #returns just first few rows of data
str(data_raw) #describes the types of data in you rmatrix


###########################
#NETWORK ANALYSES##########
###########################

min.occu<-0.4
data.filt<-data[, apply(data==0, 2, sum) < nrow(data)-min.occu*nrow(data)]

range.0to1 <- function(x){(x-min(x))/(max(x)-min(x))}
data.range <-apply(data.filt,2,range.0to1)

corr<-rcorr(data.range, type="spearman")

corr$r[upper.tri(corr$r, diag=T)] <- NA
corr$P[upper.tri(corr$P, diag=T)] <- NA 
results<-as.data.frame(cbind(which(!is.na(corr$r), arr.ind = T, useNames=F), na.omit(as.vector(corr$r)), na.omit(as.vector(corr$P))))
results[,1] <- rownames(corr$r)[results[,1]]
results[,2] <- rownames(corr$r)[results[,2]]
colnames(results) <- c("x", "y", "r", "rawP")

fdrP<-p.adjust(results$rawP, method="fdr")

results<-cbind(results, fdrP)

rcor <- 0.6
pvalue <- 0.001

results.filt <-subset(results, abs(results$r)>rcor & results$fdrP<pvalue)
results.net <- results.filt[,c(1,2,3,5)]

#Transform to a graph object
graph.transform.weights <- function (X) {
  require(igraph)
  data.tmp <- matrix(0, nrow=nrow(X), ncol=2)
  dimnames(data.tmp)[[2]] <- c("i", "j")
  data.tmp[,1] <- X[,1]
  data.tmp[,2] <- X[,2]
  data <- data.frame(data.tmp)
  graph.data <- graph.data.frame(data, directed=F)
  E(graph.data)$weight <- abs(as.numeric(X[,3]))
  E(graph.data)$rho<- as.numeric(X[,3])
  summary(graph.data)
  cat("Average degree:",ecount(graph.data)/vcount(graph.data)*2)
  return(graph.data)
}

