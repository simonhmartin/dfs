
source("DFS.R")

################################################################################
##############################  arabidopsis ######################################
################################################################################

FS <- read.table("empirical_data/Arabidopsis/arn_lyr_72.DP5MIN58MAC2.lyrata2_lyrata4_arenosa4.sfs")


### get Dfs

dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(14,14,4))

dfs2D_data <- get.DFS2D(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(14,14,4))

### plot

png("images/DFS_arabidopsis.lyr2_lyr4_arn4.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()

################################################################################
##############################  datepalms ######################################
################################################################################

FS <- read.table("empirical_data/datepalms/Flowers.SNPs.DP8.dactylifera_ME_dactylifera_NA_theophrasti.subsample10.sfs")

### get Dfs

dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(20,20,10))
D <- get.D.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(20,20,10))
# f <- get.f.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(20,20,10))
# dcfs <- get.dcfs(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(12,12,4))

### plot

png("images/DFS_Datepalms_dacME_dacNA_the.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()

################################################################################
################################  sparrows  ####################################
################################################################################

FS <- read.table("empirical_data/sparrows/Walsh2018_Evolution.RawSNPData.DP4.NEL_allo_NEL_sym_SAL_allo.subsample12.sfs")

### get Dfs

dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(12,12,4))
# D <- get.D.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(12,12,4))
# f <- get.f.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(12,12,4))
# dcfs <- get.dcfs(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(12,12,4))

### plot

png("images/DFS_sparrows_NELallo_NELsym_SALallo.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", ylim=c(-0.25,0.25), col_D="red", no_xlab=T)
dev.off()

################################################################################
###############################  heliconius  ###################################
################################################################################

pop_combos <- c("mpg_ama_txn_slvr   ",
                "chi_txn_ama_slv",
                "mpg_ros_chi_slv",
                "flo_chi_ros_slv")


#in this case the SFS was not polarised when it was computed.
#we therefore first polarise each SFS using the oygroup (slv) before computing DFS

for (pops in pop_combos){
    FS_heli <- read.table(paste0("empirical_data/Heliconius/bar92.DP8MIN92BIHET75.minor.",pops,".sfs"))

    #polarize counts as these are minor allele counts
    base_counts_heli <- polarize.counts(FS_heli[,-5], Ns=c(20,20,20,4))[,-4]

    N_heli = c(20,20,20)

    ### get Dfs

    dfs_data <- get.DFS(base_counts=base_counts_heli, site_counts=FS_heli[,5], N_heli)
    # D <- get.D.from.base.counts(base_counts=base_counts_heli, site_counts=FS_heli[,5], N_heli)
    # f <- get.f.from.base.counts(base_counts=base_counts_heli, site_counts=FS_heli[,5], N_heli)
    # dcfs <- get.dcfs(base_counts=base_counts_heli, site_counts=FS_heli[,5], N_heli)

    ### plot

    png(paste0("images/DFS_Heliconius_",pops,".png"), width = 2000, height = 1000, res=300)
    par(mar=c(1,4,1,1))
    plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
    dev.off()
    }


################################################################################
##############################  sticklebacks  ##################################
################################################################################

# FS <- read.table("empirical_data/sticklebacks/groves2019.DP5.pun_sin_tym.subsample10.sfs")
FS <- read.table("empirical_data/sticklebacks/groves2019.exclchr12.DP5.pun_sin_tym.subsample10.sfs")

### get Dfs

dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=site_counts <- FS[,4], c(10,10,10))
# D <- get.D.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(10,10,10))
# f <- get.f.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = c(10,10,10))

### plot


png("images/DFS_Stickelback_pun_sin_tym.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()


################################################################################
#################################  humans  #####################################
################################################################################


FS <- read.table("empirical_data/humans/chimp_1000G_DNK_ALT.chr1.GWD_GBR_nea.sample100.sfs")

dfs_data <- get.DFS(base_counts=FS[,-4], site_counts=site_counts <- FS[,4], c(100,100,2))
# D <- get.D.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = Ns)
# f <- get.f.from.base.counts(base_counts=FS[,-4], site_counts=FS[,4], Ns = Ns)
# dcfs <- get.dcfs(base_counts=FS[,-4], site_counts=FS[,4], Ns =Ns)

png("images/chimp_1000G_DNK_ALT.chr1.GWD_GBR_nea.sample100.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()

