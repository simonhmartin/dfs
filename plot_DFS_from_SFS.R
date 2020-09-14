
# Simon H. Martin 2020
# simon.martin@ed.ac.uk

################################ Start Here ####################################

# This script accompanies the paper:
# "Signatures of introgression across the allelel frequency spectrum"
# by Simon H. Martin and William Amos

# Each section below computes and plots the D frequency spectrum (DFS)
# from a different empirical site frequency spectrum (SFS).

# The input SFS is provided in tabular format:
# A tab-delimited table in which the first three columns (or 4 in the case of a 4D SFS)
# give the allele count in each population (equivalent to the indices of the multidimensional SFS).
# The final column gives the corresponding number of sites.

# In most cases the input SFS is 3D, with the first three columns corresponding to
# populations P1, P2 and P3 (see the paper for details). This means that the SFS is
# polarized, and the outgroup is assumed to carry the ancestral allele at these sites.

# In the case of Helcionius butterflies, we provide an example of using a 4D SFS.
# In this case, the input SFS is not polarzed. This means that the frequencies provided
# correspond to minor allele counts, and the fourth column gives the count in teh outgroup.
# Before DFS can be computed from a 4D SFS, it must first be polarized, and sites at
# which the outgroup is polymorpic ust be discarded.

#The functions to compute DFS and related statistics are provided in the accompanything script DFS.R

#First import these functions
source("DFS.R")

################################################################################
##############################  Arabidopsis ####################################
################################################################################

### import the frequency spectrum
FS <- read.table("empirical_data/Arabidopsis/arn_lyr_72.DP5MIN58MAC2.lyrata2_lyrata4_arenosa4.sfs")

### get Dfs

dfs_data <- get.DFS(base_counts=FS[,-4], #base counts are the first three columns (i.e everything minus column 4)
                    site_counts=FS[,4], # site counts are column 4
                    Ns = c(14,14,4)) # Ns provide the haploid sample sizes of each population (1 and 2 must always be equal) 

### plot

png("images/DFS_arabidopsis.lyr2_lyr4_arn4.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()

### code for exporting a table of plotted values
# write.table(data.frame(D=round(dfs_data$DFS,4),
#                        weight=round(dfs_data$weights,4)),
#             file="empirical_data/Arabidopsis/DFS_arabidopsis.lyr2_lyr4_arn4.csv",
#             sep=",",row.names=FALSE,quote=FALSE)


################################################################################
############################### Datepalms ######################################
################################################################################

### import the frequency spectrum
FS <- read.table("empirical_data/datepalms/Flowers.SNPs.DP8.dactylifera_Iraq_dactylifera_Morocco_theophrasti.subsample10.sfs")

### get Dfs

dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(20,20,10))

### plot

png("images/DFS_Datepalms_dacIrq_dacMor_the.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()

# write.table(data.frame(D=round(dfs_data$DFS,4),
#                        weight=round(dfs_data$weights,4)),
#             file="empirical_data/datepalms/DFS_Datepalms_dacIrq_dacMor_the.csv",
#             sep=",",row.names=FALSE,quote=FALSE)

################################################################################
################################  sparrows  ####################################
################################################################################

### import the frequency spectrum
FS <- read.table("empirical_data/sparrows/Walsh2018_Evolution.RawSNPData.DP4.NEL_allo_NEL_sym_SAL_allo.subsample12.sfs")

### get Dfs

dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(12,12,4))

### plot
png("images/DFS_sparrows_NELallo_NELsym_SALallo.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", ylim=c(-0.25,0.25), col_D="red", no_xlab=T)
dev.off()

# write.table(data.frame(D=round(dfs_data$DFS,4),
#                        weight=round(dfs_data$weights,4)),
#             file="empirical_data/sparrows/DFS_sparrows_NELallo_NELsym_SALallo.csv",
#             sep=",",row.names=FALSE,quote=FALSE)

################################################################################
###############################  heliconius  ###################################
################################################################################

#Here we have multiple inputs with different populations, so we first defien the population sets
pop_combos <- c("mpg_ama_txn_slv",
                "chi_txn_ama_slv",
                "mpg_ros_chi_slv",
                "flo_chi_ros_slv")

#for each combo, we compute dfs and plot
for (pops in pop_combos){
    #import frequenxy spectrum
    FS_heli <- read.table(paste0("empirical_data/Heliconius/bar92.DP8MIN92BIHET75.minor.",pops,".sfs"))
    
    #in this case the SFS was not polarised when it was computed (i.e it represents minor allele count).
    #we therefore first polarise each SFS using the ougroup (slv, column 5) before computing DFS
    base_counts_heli <- polarize.counts(FS_heli[,-5], Ns=c(20,20,20,4))[,-4]

    ### get Dfs
    
    dfs_data <- get.DFS(base_counts=base_counts_heli, site_counts=FS_heli[,5], Ns=c(20,20,20))
    
    ### plot
    
    png(paste0("images/DFS_Heliconius_",pops,".png"), width = 2000, height = 1000, res=300)
    par(mar=c(1,4,1,1))
    plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
    dev.off()
    
    ### code to export a table of plotted values
#     write.table(data.frame(D=round(dfs_data$DFS,4),
#                            weight=round(dfs_data$weights,4)),
#                 file=paste0("empirical_data/Heliconius/DFS_Heliconius_",pops,".csv"),
#                 sep=",",row.names=FALSE,quote=FALSE)
    }


################################################################################
##############################  sticklebacks  ##################################
################################################################################

### import the frequency spectrum
FS <- read.table("empirical_data/sticklebacks/groves2019.exclchr12.DP5.pun_sin_tym.subsample10.sfs")

### get Dfs
dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=site_counts <- FS[,4], c(10,10,10))

### plot
png("images/DFS_Stickelback_pun_sin_tym.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()

### code to export a table of plotted values
# write.table(data.frame(D=round(dfs_data$DFS,4),
#                        weight=round(dfs_data$weights,4)),
#             file="empirical_data/sticklebacks/DFS_Stickelback_pun_sin_tym.csv",
#             sep=",",row.names=FALSE,quote=FALSE)


################################################################################
#################################  humans  #####################################
################################################################################


### import the frequency spectrum
FS <- read.table("empirical_data/humans/chimp_1000G_DNK_ALT.chr1.GWD_GBR_nea.sample100.sfs")

### get Dfs
dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=site_counts <- FS[,4], Ns = c(100,100,2))

# We can also compute related statistics
D <- get.D.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(100,100,2)) #overall D
f <- get.f.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(100,100,2)) # overall fraction of introgression
dcfs <- get.dcfs(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(100,100,2)) # dcfs from Yang et al. 2012

### plot
png("images/chimp_1000G_DNK_ALT.chr1.GWD_GBR_nea.sample100.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()

### code to export a table of plotted values
# write.table(data.frame(D=round(dfs_data$DFS,4),
#                        weight=round(dfs_data$weights,4)),
#             file="empirical_data/humans/DFS_chimp_1000G_DNK_ALT.chr1.GWD_GBR_nea.sample100.csv",
#             sep=",",row.names=FALSE,quote=FALSE)
