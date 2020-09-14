
get.patterns <- function(p1, p2, p3, pO=NULL){
    #assume that if there's no outgroup, it is always fixed ancestral
    if (is.null(pO)==TRUE | length(pO) == 0) pO <- 0
    
    data.frame(ABBA = (1 - p1) * p2 * p3 * (1-pO),
               BABA = p1 * (1 - p2) * p3 * (1-pO),
               ABBA_BAAB = (1 - p1) * p2 * p3 * (1-pO) + p1 * (1-p2) * (1-p3) * pO,
               BABA_ABAB = p1 * (1 - p2) * p3 * (1-pO) + (1-p1) * p2 * (1-p3) * pO,
               ABBA_f = (1 - p1) * p3 * p3 * (1-pO),
               BABA_f = p1 * (1-p3) * p3 * (1-pO),
               ABBA_BAAB_f = (1 - p1) * p3 * p3 * (1-pO) + p1 * (1-p3)**2 * pO,
               BABA_ABAB_f = p1 * (1 - p3) * p3 * (1-pO) + (1-p1) * p3 * (1-p3) * pO
               )
    }



get.DFS <- function(base_counts, site_counts=NULL, Ns=NULL){
    
    #if the number of haplotypes per population is not specified, assume it is the maximum value observed - NOT RECOMMENDED
    if (is.null(Ns) == TRUE) Ns <- apply(base_counts, 2, max)
    
    #site counts are used to compress the input data. They give the number of sites observed with those base counts
    #if not specified, assume all patterns are represented once
    if (is.null(site_counts) == TRUE) site_counts <- rep(1, nrow(base_counts))
    
    if (!(((ncol(base_counts) == 4 & length(Ns) ==4) | (ncol(base_counts) == 3 & length(Ns) == 3)) & Ns[1] == Ns[2])){
        print("Incorrect specification")
        return()
    }
    
    #convert base_counts to frequencies
    freqs = base_counts / t(replicate(nrow(base_counts), Ns))
    
    N = Ns[1]
    
    #identify sites where P1 and P2 each have a specified allele frequency
    idx_1 = lapply(1:N, function(i) which(base_counts[,1] == i))
    idx_2 = lapply(1:N, function(i) which(base_counts[,2] == i))
    
    #get total counts of each pattern
    pattern_sums_by_count_1 <- sapply(1:N, function(i) apply(get.patterns(freqs[idx_1[[i]],1],
                                                                        freqs[idx_1[[i]],2],
                                                                        freqs[idx_1[[i]],3],
                                                                        freqs[idx_1[[i]],4]) * site_counts[idx_1[[i]]],2,sum))

    pattern_sums_by_count_2 <- sapply(1:N, function(i) apply(get.patterns(freqs[idx_2[[i]],1],
                                                                        freqs[idx_2[[i]],2],
                                                                        freqs[idx_2[[i]],3],
                                                                        freqs[idx_2[[i]],4]) * site_counts[idx_2[[i]]],2,sum))
    
    ABBA_by_count <- pattern_sums_by_count_2["ABBA",]
    BABA_by_count <- pattern_sums_by_count_1["BABA",]
    
    DFS <- (ABBA_by_count - BABA_by_count) /
            (ABBA_by_count + BABA_by_count)

    weights <-     (ABBA_by_count + BABA_by_count)/
               sum((ABBA_by_count + BABA_by_count))
    
    data.frame(DFS=DFS, weights=weights, ABBA=ABBA_by_count, BABA=BABA_by_count)
    }



get.D.from.base.counts <- function(base_counts, site_counts=NULL, Ns=NULL, full=FALSE){
    
    #if the number of haplotypes per population is not specified, assume it is the maximum value observed - NOT RECOMMENDED
    if (is.null(Ns) == TRUE) Ns <- apply(base_counts, 2, max)
    
    #site counts are used to compress the input data. They give the number of sites observed with those base counts
    #if not specified, assume all patterns are represented once
    if (is.null(site_counts) == TRUE) site_counts <- rep(1, nrow(base_counts))
    
    if (!((ncol(base_counts) == 4 & length(Ns) ==4) | (ncol(base_counts) == 3 & length(Ns) == 3)) ){
        print("Incorrect specification")
        return()
    }
    
    #convert base_counts to frequencies
    freqs = base_counts / t(replicate(nrow(base_counts), Ns))

    idx <- 1:nrow(freqs)
    
    #get total counts of each pattern
    pattern_sums <- apply(get.patterns(freqs[idx,1],freqs[idx,2],freqs[idx,3],freqs[idx,4]) * site_counts,2,sum)
        
    if (full == FALSE){
        ABBA <- pattern_sums["ABBA"]
        BABA <- pattern_sums["BABA"]
        return(as.numeric((ABBA - BABA) / (ABBA + BABA)))
        }
    else{
        ABBA_BAAB <- pattern_sums["ABBA_BAAB"]
        BABA_ABAB <- pattern_sums["BABA_ABAB"]
        return(as.numeric((ABBA_BAAB - BABA_ABAB) / (ABBA_BAAB + BABA_ABAB)))
        }
    }

get.f.from.base.counts <- function(base_counts, site_counts=NULL, Ns=NULL, full=FALSE){
    
    #if the number of haplotypes per population is not specified, assume it is the maximum value observed - NOT RECOMMENDED
    if (is.null(Ns) == TRUE) {
        print("WARNING: Ns not provided, assuming maximum count observed per population is N haploid samples")
        Ns <- apply(base_counts, 2, max)
        }
    
    #site counts are used to compress the input data. They give the number of sites observed with those base counts
    #if not specified, assume all patterns are represented once
    if (is.null(site_counts) == TRUE) site_counts <- rep(1, nrow(base_counts))
    
    if (!((ncol(base_counts) == 4 & length(Ns) ==4) | (ncol(base_counts) == 3 & length(Ns) == 3)) ){
        print("Incorrect specification")
        return()
    }
    
    #convert base_counts to frequencies
    freqs = base_counts / t(replicate(nrow(base_counts), Ns))

    idx <- 1:nrow(freqs)
    
    #get total counts of each pattern
    pattern_sums <- apply(get.patterns(freqs[idx,1],freqs[idx,2],freqs[idx,3],freqs[idx,4]) * site_counts,2,sum)
    
    if (full == FALSE){
        ABBA <- pattern_sums["ABBA"]
        BABA <- pattern_sums["BABA"]
        ABBA_f <- pattern_sums["ABBA_f"]
        BABA_f <- pattern_sums["BABA_f"]
        return(as.numeric((ABBA - BABA) / (ABBA_f - BABA_f)))
        }
    else{
        ABBA_BAAB <- pattern_sums["ABBA_BAAB"]
        BABA_ABAB <- pattern_sums["BABA_ABAB"]
        ABBA_BAAB_f <- pattern_sums["ABBA_BAAB_f"]
        BABA_ABAB_f <- pattern_sums["BABA_ABAB_f"]
        return(as.numeric((ABBA_BAAB - BABA_ABAB) / (ABBA_BAAB_f - BABA_ABAB_f)))
        }
    }


get.f4.from.base.counts <- function(base_counts, site_counts=NULL, Ns=NULL) {
    #if the number of haplotypes per population is not specified, assume it is the maximum value observed - NOT RECOMMENDED
    if (is.null(Ns) == TRUE) Ns <- apply(base_counts, 2, max)
    
    freqs = base_counts / t(replicate(nrow(base_counts), Ns))
    f4_by_site <- (freqs[,1]-freqs[,2])*(freqs[,3]-freqs[,4])
    if(is.null(site_counts) == TRUE) return(mean(f4_by_site))
    else weighted.mean(f4_by_site,site_counts)
    }


get.dcfs <- function(base_counts, site_counts=NULL, Ns=NULL){
    
    #if the number of haplotypes per population is not specified, assume it is the maximum value observed - NOT RECOMMENDED
    if (is.null(Ns) == TRUE) Ns <- apply(base_counts, 2, max)
    
    #site counts are used to compress the input data. They give the number of sites observed with those base counts
    #if not specified, assume all patterns are represented once
    if (is.null(site_counts) == TRUE) site_counts <- rep(1, nrow(base_counts))
    
    #assert that base counts are for three populations (they must be polarized)
    if (ncol(base_counts) != 3 | length(Ns) != 3){
        print("Incorrect specification")
        return()
    }
    
    #convert base_counts to frequencies
    freqs = base_counts / t(replicate(nrow(base_counts), Ns))
    
    N = Ns[2]
    
    #identify sites P2 has a specified allele frequency
    indices = lapply(1:N, function(i) which(base_counts[,2] == i))
    
    #for each frequency in P2, multiply by probability of getting an ancestral allele in P1 and derived in P3, and sum site counts
    dcfs_unscaled <- sapply(indices, function(idx) sum((1-freqs[idx,1]) * freqs[idx,3] * site_counts[idx])) 
    
    dcfs_unscaled/sum(dcfs_unscaled)
    }

get.DFS2D <- function(base_counts, site_counts=NULL, Ns=NULL){
    
    #if the number of haplotypes per population is not specified, assume it is the maximum value observed - NOT RECOMMENDED
    if (is.null(Ns) == TRUE) Ns <- apply(base_counts, 2, max)
    
    #site counts are used to compress the input data. They givethe number of sites observed with those base counts
    #if not specified, assume all patterns are represented once
    if (is.null(site_counts) == TRUE) site_counts <- rep(1, nrow(base_counts))
    
    if (!(((ncol(base_counts) == 4 & length(Ns) ==4) | (ncol(base_counts) == 3 & length(Ns) == 3)) & Ns[1] == Ns[2])){
        print("Incorrect specification")
        return()
    }
    
    #convert base_counts to frequencies
    freqs = base_counts / t(replicate(nrow(base_counts), Ns))
    
    #identify sites where P1 and P2 each have a specified allele frequency
    pattern_sums_array <- array(dim=c(Ns[1],Ns[2],2), dimnames = list(1:Ns[1], 1:Ns[2], c("ABBA", "BABA")))
    
    for (i in 1:Ns[1]){
        for (j in 1:Ns[2]){
            idx = which(base_counts[,1] == i & base_counts[,2] == j)
            pattern_sums_array[i,j,] <- apply(get.patterns(freqs[idx,1],
                                                           freqs[idx,2],
                                                           freqs[idx,3],
                                                           freqs[idx,4])[,c("ABBA","BABA")] * site_counts[idx],2,sum)
            }
        }
    
    DFS <- (pattern_sums_array[,,"ABBA"] - pattern_sums_array[,,"BABA"]) /
            (pattern_sums_array[,,"ABBA"] + pattern_sums_array[,,"BABA"])

    weights <- (pattern_sums_array[,,"ABBA"] + pattern_sums_array[,,"BABA"])/
            sum((pattern_sums_array[,,"ABBA"] + pattern_sums_array[,,"BABA"]))
    
    list(DFS2D=DFS, weights=weights)
    }



polarize.counts <- function(counts, Ns, OGcolumn=NULL, outgroup_pol_to_NA=TRUE){
    #If outgroup column is not specified, assume it is last one 
    OG <- ifelse(is.null(OGcolumn) == TRUE, ncol(counts), OGcolumn)
    counts_plr <- counts
    if (outgroup_pol_to_NA == TRUE){
        counts_plr[counts_plr[,OG]!=0 & counts_plr[,OG]!=Ns[OG],] <- NA
        }
    #flip sites where numbers must be flipped
    flip_idx = which(counts_plr[,OG] > Ns[OG]/2)
    counts_plr[flip_idx,] <- t(apply(counts[flip_idx,], 1, function(x) Ns-x))
    counts_plr
    }


################################################ plotting functions

plotDFS <- function(DFS, weights, method="lines", ylim=c(-1,1), show_D=TRUE,
                    col="black", col_D="black", width_scale=100, no_xlab=FALSE, add=FALSE){
    
    if (method == "lines"){
        N = length(DFS)
        if (add == FALSE){
            plot(0, xlim = c(1,N), ylim = ylim, cex=0, xlab = "", ylab = "", xaxt="n", bty="n")
            abline(h=0)
            }
        segments(1:N, 0, 1:N, DFS, lwd = width_scale*weights, lend=1, col=col)
        }
    
    if (method == "bars") barplot(DFS, col= rgb(0,0,0,weights), ylim = ylim, add=add)
    
    if (method == "scaled_bars") barplot(DFS*weights, ylim = ylim, add=add)
    
    if (no_xlab == FALSE & add == FALSE) mtext(1,text="Derived allele frequency", line = 0)
    
    if (add == FALSE) mtext(2,text=expression(italic("D")), line = 2.8, las=2)
    if (show_D == TRUE) abline(h= sum(DFS * weights), lty = 2, col=col_D)
    
    }


plot.dcfs <- function(dcfs){
    plot(dcfs, type="b")
    }


################################################ other miscilaneous functions

#Most function here use a table format for the SFS, this is more suitable for a large, sparce array
#This function converts to the more conventional NxN(xN...) array.

sfs.table.to.array <- function(sfs_table, dims=NULL, count_col=NULL){
    if (is.null(dims)==TRUE) dims <- apply(sfs_table[,-ncol(sfs_table)], 2, max) + 1
    ndim <- length(dims)
    if (is.null(count_col)==TRUE) count_col <- ndim+1
    arr <- array(dim=dims)
    for (i in 1:nrow(sfs_table)){
        arr[as.matrix(sfs_table[i,-count_col])+1] <- sfs_table[i,count_col]
        }
    arr
    }

