
# Simon H. Martin 2020
# simon.martin@ed.ac.uk

################################ Start Here ####################################

# This script accompanies the paper:
# "Signatures of introgression across the allelel frequency spectrum"
# by Simon H. Martin and William Amos

# It provides code for setting up and simulating site frequency spectra
# under various models, and then plotting the D frequency spectrum (DFS).

# Simulations are performed using moments (Jouganous et al. 2017).

#import necessary libraries
library(parallel)
library(data.table)

#import code for computaing and plotting DFS and related statistics
source("DFS.R")

#change limit for scientific notation (this avoids exporting numbers in strange fromat when making commands)
options(scipen = 10)

################################################################################
############# Functions for making the simulation input commands ###############
################################################################################


### NOTE by default in moments (as in dadi), first split produces P2 from P1. Second produced P3 from P2
### Thus, for classic ABBA BABA, P1 and P3 are reversed (e.g. P1 is Neanderthal and P3 is Africa)

if.null <- function(x,y){
    if (is.null(x)==TRUE) return(y)
    return(x)
    }

#make moments simulation command for three population models
make_3pop_sim_command <- function(params){
    paste("python", "sim_SFS_moments.py",
        "--Nsam", paste(params$n1,params$n2,params$n3),
        
        "--twoPopPeriod", paste0("T=",params$p2_T, #input is in units of 2N gen
                                 " N1=",if.null(params$p2_N1,1),
                                 " N2=",if.null(params$p2_N2,1),
                                 " m12=",if.null(params$p2_m12,0),
                                 " m21=",if.null(params$p2_m21,0)),
        
        "--threePopPeriod", paste0("T=",params$p31_T,
                                 " N1=",if.null(params$p31_N1,1),
                                 " N2=",if.null(params$p31_N2,1),
                                 " N3=",if.null(params$p31_N3,1),
                                 " m12=",if.null(params$p31_m12,0),
                                 " m21=",if.null(params$p31_m21,0),
                                 " m23=",if.null(params$p31_m23,0),
                                 " m32=",if.null(params$p31_m32,0),
                                 " m13=",if.null(params$p31_m13,0),
                                 " m31=",if.null(params$p31_m31,0)),
    
        "--threePopPeriod", paste0("T=",params$p32_T,
                                 " N1=",if.null(params$p32_N1,1),
                                 " N2=",if.null(params$p32_N2,1),
                                 " N3=",if.null(params$p32_N3,1),
                                 " m12=",if.null(params$p32_m12,0),
                                 " m21=",if.null(params$p32_m21,0),
                                 " m23=",if.null(params$p32_m23,0),
                                 " m32=",if.null(params$p32_m32,0),
                                 " m13=",if.null(params$p32_m13,0),
                                 " m31=",if.null(params$p32_m31,0)),
    
        " 2> /dev/null")
    }

#make moments simulation command for four population models
make_4pop_sim_command <- function(params){
    paste("python", "sim_SFS_moments.py",
        "--Nsam", paste(params$n1,params$n2,params$n3,params$n4),
        
        "--twoPopPeriod", paste0("T=",params$p2_T, #input is in units of 2N gen
                                 " N1=",if.null(params$p2_N1,1),
                                 " N2=",if.null(params$p2_N2,1),
                                 " m12=",if.null(params$p2_m12,0),
                                 " m21=",if.null(params$p2_m21,0)),
        
        "--threePopPeriod", paste0("T=",params$p3_T,
                                 " N1=",if.null(params$p3_N1,1),
                                 " N2=",if.null(params$p3_N2,1),
                                 " N3=",if.null(params$p3_N3,1),
                                 " m12=",if.null(params$p3_m12,0),
                                 " m21=",if.null(params$p3_m21,0),
                                 " m13=",if.null(params$p3_m13,0),
                                 " m31=",if.null(params$p3_m31,0),
                                 " m23=",if.null(params$p3_m23,0),
                                 " m32=",if.null(params$p3_m32,0)),
        
        "--fourPopPeriod", paste0("T=",params$p41_T,
                                 " N1=",if.null(params$p41_N1,1),
                                 " N2=",if.null(params$p41_N2,1),
                                 " N3=",if.null(params$p41_N3,1),
                                 " N4=",if.null(params$p41_N4,1),
                                 " m12=",if.null(params$p41_m12,0),
                                 " m21=",if.null(params$p41_m21,0),
                                 " m13=",if.null(params$p41_m13,0),
                                 " m31=",if.null(params$p41_m31,0),
                                 " m14=",if.null(params$p41_m14,0),
                                 " m41=",if.null(params$p41_m41,0),
                                 " m23=",if.null(params$p41_m23,0),
                                 " m32=",if.null(params$p41_m32,0),
                                 " m24=",if.null(params$p41_m24,0),
                                 " m42=",if.null(params$p41_m42,0),
                                 " m34=",if.null(params$p41_m34,0),
                                 " m43=",if.null(params$p41_m43,0)),
        
        "--fourPopPeriod", paste0("T=",params$p42_T,
                                 " N1=",if.null(params$p42_N1,1),
                                 " N2=",if.null(params$p42_N2,1),
                                 " N3=",if.null(params$p42_N3,1),
                                 " N4=",if.null(params$p42_N4,1),
                                 " m12=",if.null(params$p42_m12,0),
                                 " m21=",if.null(params$p42_m21,0),
                                 " m13=",if.null(params$p42_m13,0),
                                 " m31=",if.null(params$p42_m31,0),
                                 " m14=",if.null(params$p42_m14,0),
                                 " m41=",if.null(params$p42_m41,0),
                                 " m23=",if.null(params$p42_m23,0),
                                 " m32=",if.null(params$p42_m32,0),
                                 " m24=",if.null(params$p42_m24,0),
                                 " m42=",if.null(params$p42_m42,0),
                                 " m34=",if.null(params$p42_m34,0),
                                 " m43=",if.null(params$p42_m43,0)),
    
        " 2> /dev/null")
    }

#function to convert parameter names to those used by the moments script.
#This is just because, like dadi, moments numbers pops differently from the typical ABBA BABA numbering
make_moments_params_3pop <- function(params){
                  list(n1 = params$n3,
                       n2 = params$n2,
                       n3 = params$n1, #number of samples in each pop
                       
                       p2_T = params$p2_T, #length of two pop periods (in 2N generations)
                       p2_m21 = params$p2_m12_3, #first three pop period m23
                       p2_m12 = params$p2_m3_12, #first three pop period m32
                       p2_N1 = params$p2_N3, #first twp pop period N1
                       p2_N2 = params$p2_N12, #first twp pop period N1
                       
                       p31_T =   params$p31_T, #first three pop period
                       p31_N1 =  params$p31_N3, #first three pop period N1
                       p31_N2 =  params$p31_N2, #first three pop period N2
                       p31_N3 =  params$p31_N1, #first three pop period N2
                       p31_m21 = params$p31_m23, #first three pop period m23
                       p31_m12 = params$p31_m32, #first three pop period m32
                       p31_m31 = params$p31_m13, #first three pop period m31
                       p31_m13 = params$p31_m31, #first three pop period m13
                       p31_m23 = params$p31_m21, #first three pop period m23
                       p31_m32 = params$p31_m12, #first three pop period m32
                       
                       p32_T = params$p32_T, #second three pop period
                       p32_N1 = params$p32_N3, #second three pop period N1
                       p32_N2 = params$p32_N2, #second three pop period N2
                       p32_N3 = params$p32_N1, #second three pop period N2
                       p32_m21 = params$p32_m23, #second three pop period m23
                       p32_m12 = params$p32_m32, #second three pop period m32
                       p32_m31 = params$p32_m13, #second three pop period m31
                       p32_m13 = params$p32_m31, #second three pop period m13
                       p32_m23 = params$p32_m21, #second three pop period m23
                       p32_m32 = params$p32_m12 #second three pop period m32
            )
        }

        
##################################################################################
##################################################################################
####################  RUN SINGLE SIMULATION AND PLOT  ############################
##################################################################################
##################################################################################

#below is code for running a single simulation and plotting the result directly

##################################################################################
#########################  standard 3 pop model  #################################
##################################################################################

# NOTE in moments, as in in dadi, first split produces pop 2 from 1. Second produces 3 from 2
# So for ABBA BABA traditional P3 (e.g. Neanderthal) is pop 1 in dadi and P1 is pop 3 in moments
# so the migration is simulated as between pops 1 and 2

# This was dealt with above using the make_moments_params_3pop function,
# but here we just remember toswitch the numbers

#Three population model
params <- list(n1 = 4, n2=20, n3=20, #number of samples in each pop
               
               p2_T = 0.2, #length of two pop periods (in 2N generations)
               p2_m21 = 0, #first three pop period m23
               p2_m12 = 0, #first three pop period m32
               
               p31_T = 0.1, #first three pop period
               p31_N1 = 1, #first three pop period N1
               p31_N2 = 1, #first three pop period N2
               p31_N3 = 1, #first three pop period N2
               p31_m21 = 0, #first three pop period m23
               p31_m12 = 0, #first three pop period m32
               p31_m31 = 0, #first three pop period m31
               p31_m13 = 0, #first three pop period m13

               p32_T = 0.1, #second three pop period
               p32_N1 = 1, #second three pop period N1
               p32_N2 = 1, #second three pop period N2
               p32_N3 = 1, #second three pop period N3
               p32_m21 = 1, #second three pop period m23
               p32_m12 = 0, #second three pop period m32
               p32_m31 = 0, #second three pop period m31
               p32_m13 = 0, #second three pop period m13
               p32_m23 = 0, #second three pop period m31
               p32_m32 = 0 #second three pop period m13
               )

#population indices (remember dadi first split produced P2 from P1 and second split produced P3 from P2)
P1 <- 3
P2 <- 2
P3 <- 1

#make command
command <- make_3pop_sim_command(params)

print(command)

#run simulation
SFS <- read.table(text=system(command, intern=TRUE,ignore.stderr=FALSE))

#get DFS
dfs_data <- get.DFS(base_counts=SFS[,c(P1,P2,P3)], site_counts=SFS[,ncol(SFS)], Ns=c(params$n1,params$n2,params$n3,params$n4)[c(P1,P2,P3)])

#plot
plotDFS(dfs_data$DFS, dfs_data$weights, col_D="red")

dfs_data

##################################################################################
#######################  ancestral structure model  ##############################
##################################################################################


# NOTE in moments, first split produces dadi pop 2 from 1. Second produces 3 from 2
# for ancestral structure model, this is what we want
# because the first split introduces structure (e.g. in Africa) and second splits out pop off (e.g. neanderthal)  (see Figure in Yang et al. 2012)

params <- list(n1 = 20, n2=20, n3=4, #number of samples in each pop
               
               p2_T = 1, #length of ancestral structure period pop period (in N generations)
               p2_m12 = 2, #migration during ancestral structure period
               p2_m21 = 2, #migration during ancestral structure period
               
               p31_T = 0.1, #length of continued structure after first split
               p31_m12 = 2, #first three pop period m23
               p31_m21 = 2, #first three pop period m23
               
               p32_T = 0.1 #split between ingroups
               )

#population indices for this model the pops are correct because first split produces e.g. Neanderthal, which is P1 here
P1 <- 1
P2 <- 2
P3 <- 3

#make command
command <- make_3pop_sim_command(params)

print(command)

#run simulation
SFS <- read.table(text=system(command, intern=TRUE,ignore.stderr=FALSE))


#get DFS
dfs_data <- get.DFS(base_counts=SFS[,c(P1,P2,P3)], site_counts=SFS[,ncol(SFS)], Ns=c(params$n1,params$n2,params$n3,params$n4)[c(P1,P2,P3)])

# for comparison, also get dcfs (Yang et al. 2012)
dcfs <- get.dcfs(base_counts=SFS[,c(P1,P2,P3)], site_counts=SFS[,4], Ns=c(params$n1,params$n2,params$n3)[c(P1,P2,P3)])

#plot both DFS and dcfs

par(mfrow = c(2,1))

plotDFS(dfs_data$DFS, dfs_data$weights)

plot(dcfs, type="b")


##################################################################################
#################################  4 pop model  ##################################
##################################################################################

params <- list(n1 = 4, n2=10, n3=10, n4=10, #number of samples in each pop
               
               p2_T = 0.1, #length of two pop periods (in 2N generations)
               
               p3_T = 0.1, #first three pop period
               p3_N1 = 1, #first three pop period N1
               p3_N2 = 1, #first three pop period N2
               p3_N3 = 1, #first three pop period N2
               p3_m21 = 0, #first three pop period m23
               p3_m12 = 0, #first three pop period m32
               p3_m31 = 0, #first three pop period m31
               p3_m13 = 0, #first three pop period m13
               
               p41_T = 0.1, #second three pop period
               p41_N1 = 1, #second three pop period N1
               p41_N2 = 1, #second three pop period N2
               p41_N3 = 1, #second three pop period N3
               p41_N4 = 1, #second three pop period N4
               p41_m23 = 0, #second three pop period m31
               p41_m32 = 0, #second three pop period m13
               p41_m24 = 0, #second three pop period m31
               p41_m42 = 0, #second three pop period m13
               p41_m34 = 0, #second three pop period m31
               p41_m43 = 0, #second three pop period m13
               
               p42_T = 0.1, #second three pop period
               p42_N1 = 1, #second three pop period N1
               p42_N2 = 1, #second three pop period N2
               p42_N3 = 1, #second three pop period N3
               p42_N4 = 1, #second three pop period N4
               p42_m23 = 0, #second three pop period m31
               p42_m32 = 1, #second three pop period m13
               p42_m24 = 0, #second three pop period m31
               p42_m42 = 0, #second three pop period m13
               p42_m34 = 0, #second three pop period m31
               p42_m43 = 0 #second three pop period m13
               )

#population indices (remember in moments first split produced P2 from P1 and second split produced P3 from P2)
P1 <- 4
P2 <- 3
P3 <- 2
OG <- 1

command <- make_4pop_sim_command(params)

#run simulation
SFS <- read.table(text=system(command, intern=TRUE,ignore.stderr=FALSE))


# With four populations there are different ways to compute DFS.
# In this case, we have a polarised SFS, so we can add the fourth population and compute as before.
# This way, we will also be using sites that are polymorphic in the outgroup, but scaling according to the frequency of thenacestral allele

dfs_data_OG <- get.DFS(base_counts=SFS[,c(P1,P2,P3,OG)], site_counts=SFS[,5], Ns=c(params$n1,params$n2,params$n3,params$n4)[c(P1,P2,P3,OG)])
plotDFS(dfs_data_OG$DFS, dfs_data$weights)


# Alternatively, we can pretend that we don't know what the ancestral allele is (ie. an unpolarised SFS)
# In this case, we have to first polarise the counts, using the outgroup as a reference for what is ancestral.
# Sites that are polymorphic in the outgroup will be ignored as they are ambiguous

dfs_data_plr <- get.DFS(base_counts=polarize.counts(SFS[,-5],
                                                    Ns=c(params$n1,params$n2,params$n3,params$n4),
                                                    OGcolumn=OG,
                                                    outgroup_pol_to_NA=TRUE)[,c(P1,P2,P3,OG)],
                        site_counts=SFS[,5],
                        Ns=c(params$n1,params$n2,params$n3,params$n4)[c(P1,P2,P3,OG)])

plotDFS(dfs_data_plr$DFS, dfs_data_plr$weights)


##################################################################################
##################################################################################
############################    BATCH SIMULATIONS     ############################
##################################################################################
##################################################################################

# Below is code for simulating many thousands of combinations of parameters to explore parameter space.

expand.grid.restricted <- function(defaults, groups, max_alt_per_group, alternatives){
    #this function is like expand.grid, except if you have too many combos, you can restrict the number that can vary at any time.
    #you define groups of variables and set the number that can vary per group
    #first check that the number of groups provided matches the length of the value list
    if (length(defaults) != length(groups)){
        print("Mismatched lengths of inputs.")
        }
    #get unique list of groups
    GROUPS = unique(groups)
    if (length(GROUPS) != length(max_alt_per_group)){
        print("'max_alt_per_group' must have one value for each group.")
        return()
        }
    if (is.null(names(defaults))) names(defaults) <- 1:length(defaults)
    value_names <- names(defaults)
    #make lists of groups to vary at any given time
    groups_to_vary <- GROUPS[which(max_alt_per_group >= 1)]
    group_combos_to_vary <- list()
    for (x in 1:length(groups_to_vary)) group_combos_to_vary <- c(group_combos_to_vary, combn(groups_to_vary, x, simplify=F))
    
    #make lists of which variables will vary at any given time
    value_name_combos_to_vary_by_group <- list()
    for (group in GROUPS){
        value_name_combos_to_vary_by_group[[group]] <- list()
        if (max_alt_per_group[group] >= 1){
            for (x in 1:max_alt_per_group[group]){
                value_name_combos_to_vary_by_group[[group]] <- c(value_name_combos_to_vary_by_group[[group]],
                                                           combn(value_names[which(groups==group)], x, simplify=F))
                }
            }
        }
    #now we know which groups will vary, and which values will vary within each group
    #now make list of values to expand for each of these case
    value_lists_to_expand <- list(as.list(defaults)) #the first is all defaults only
    n=1
    for (group_combo in group_combos_to_vary){
        #now we know which groups will vary, we need to pull out the names of the values to vary from each group. There are multiple combinations of these
        group_combo_value_combos <- expand.grid(lapply(value_name_combos_to_vary_by_group, function(x) 1:length(x))[group_combo])
        #and make a vectore of the value names to vary
        for (i in 1:nrow(group_combo_value_combos)){
            current_value_names_to_vary <- unlist(lapply(1:length(group_combo),
                                                    function(j) value_name_combos_to_vary_by_group[[group_combo[j]]][[group_combo_value_combos[i,j]]]))
            #now make a list of defaults and then change those that need to be expanded
            current_value_list <- as.list(defaults)
            for (value_name in current_value_names_to_vary){
                current_value_list[[value_name]] <- alternatives[[value_name]]
                }
            n <- n+1
            value_lists_to_expand[[n]] <- current_value_list
            }
        }
    as.data.frame(rbindlist(lapply(value_lists_to_expand, expand.grid)))
    }


##################################################################################
################    set parameters (user to make changes here)   #################
##################################################################################

#set parameters common to all models

### NOTE in dadi, first split produces P2 from P1. Second produced P3 from P2
### Thus, for classic ABBA BABA, P1 and P3 are reversed (P1 is Neanderthal and P3 is Africa)

defaults <- list(n1 = 20, n2=20, n3=4, #number of samples in each pop
               
               p2_T = 0.1, #length of two pop periods (in N generations)
               p31_T = 0.1, #first three pop period
               p32_T = 0.1, #second three pop period
               
               p2_N12 = 1, #pop size ancestor of P1 and P2
               p2_N3 = 1, #pop size ancestor P3
               p2_m12_3 = 0, #migration
               p2_m3_12 = 0, #migration
                
               p31_N1 = 1, #first three pop period N1
               p31_N2 = 1, #first three pop period N2
               p31_N3 = 1, #first three pop period N2
                
               p32_N1 = 1, #second three pop period N1
               p32_N2 = 1, #second three pop period N2
               p32_N3 = 1, #second three pop period N3
               
               p31_m12 = 0, #first three pop period m12
               p31_m21 = 0, #first three pop period m21
               p31_m23 = 0, #first three pop period m23
               p31_m32 = 0, #first three pop period m32
               p31_m13 = 0, #first three pop period m13
               p31_m31 = 0, #first three pop period m31
                
               p32_m12 = 0, #second three pop period m12
               p32_m21 = 0, #second three pop period m21
               p32_m23 = 0, #second three pop period m23
               p32_m32 = 0, #second three pop period m32
               p32_m13 = 0, #second three pop period m13
               p32_m31 = 0  #second three pop period m31
               )


#set alternatives for split times, population sizes and migration rates
alt_times <- c(0.2,0.5)
alt_Ns <- c(0.1,5)
alt_migs <- c(0.5,1,2,5)

alternatives <- list(n1 = NA, n2=NA, n3=NA, #number of samples in each pop
               
               p2_T = alt_times, #length of two pop periods (in N generations)
               p31_T = alt_times, #first three pop period
               p32_T = alt_times, #second three pop period
               
               p2_N12 = alt_Ns, #pop size ancestor of P1 and P2
               p2_N3 = alt_Ns, #pop size ancestor P3
                
               p31_N1 = alt_Ns, #first three pop period N1
               p31_N2 = alt_Ns, #first three pop period N2
               p31_N3 = alt_Ns, #first three pop period N2
                
               p32_N1 = alt_Ns, #second three pop period N1
               p32_N2 = alt_Ns, #second three pop period N2
               p32_N3 = alt_Ns, #second three pop period N3
                
               p2_m12_3 = alt_migs, #two pop period m12
               p2_m3_12 = alt_migs, #two pop period m12
               
               p31_m12 = alt_migs, #first three pop period m12
               p31_m21 = alt_migs, #first three pop period m21
               p31_m23 = alt_migs, #first three pop period m23
               p31_m32 = alt_migs, #first three pop period m32
               p31_m13 = alt_migs, #first three pop period m13
               p31_m31 = alt_migs, #first three pop period m31
                
               p32_m12 = alt_migs, #second three pop period m12
               p32_m21 = alt_migs, #second three pop period m21
               p32_m23 = alt_migs, #second three pop period m23
               p32_m32 = alt_migs, #second three pop period m32
               p32_m13 = alt_migs, #second three pop period m13
               p32_m31 = alt_migs  #second three pop period m31
               )

#use date as the name for this simulation run
params_name <- Sys.Date()

##################################################################################
#######################     prepare simulations      #############################
##################################################################################

### make parameter grids

# Here we make all combinations of parameters.
# The expand.grid.restricted command takes the defaults and alternatives (two lists of same length and order)
# and also 'groups' which tells us which parameters are allowed to vary in each run.
# For example group seven has a maximum alternatives of 2
# which means at most 2 of these parameters can vary from the defaults in a given run

params_grid <-    expand.grid.restricted(defaults=defaults,
                                        groups=c(1,2,3, #samples sizes
                                                 4,5,6, #split times
                                                 7,7,7,7,7,7,7,7, #pop sizes
                                                 7,7,7,7,7,7,7,7,7,7,7,7,7,7), #migration rates
                                        max_alt_per_group=c(0,0,0,1,1,1,2),
                                        alternatives=alternatives)

# make the simulation commands that will be fed to the python script
commands <- sapply(1:nrow(params_grid), function(i) make_3pop_sim_command(make_moments_params_3pop(params_grid[i,])))

##################################################################################
#########################  RUN SIMULATIONS  ######################################
##################################################################################

# do all reps at once. might take quite a bit of RAM, but sfs is not too lareg
freqs_text <- mclapply(commands, system, intern=TRUE,ignore.stderr=FALSE, mc.cores = 20)

#Read the text in as SFS tables
counts <- lapply(freqs_text, function(t) read.table(text=t))

##################################################################################
############################## GET DFS  ##########################################
##################################################################################


#get Dfs directly from derived allele counts
# here we switch populations 1 and 3 because in moments the order of splitting is opposite to our model
dfs_data <- mclapply(1:length(counts), function(i) get.DFS(base_counts=counts[[i]][,c(3,2,1)],
                                                         site_counts=counts[[i]][,4],
                                                         Ns=as.numeric(params_grid[i,c("n1","n2","n3")])), mc.cores=20)


##################################################################################
#########################  Export Output  ########################################
##################################################################################

output <- data.frame(param_set=unlist(lapply(1:length(dfs_data), function(i) rep(i,nrow(dfs_data[[i]])))),
                     DFS=round(unlist(lapply(dfs_data, function(df) df$DFS)),4),
                     weights=round(unlist(lapply(dfs_data, function(df) df$weights)),4))

#Dfs data all in one file
write.table(output, file = paste("sim_DFS_moments", params_name, "tsv", sep="."), quote=F, row.names=F, sep="\t")

#export parameter grid
write.table(cbind(data.frame(param_set=1:nrow(params_grid)),params_grid),
                  file = paste("params", params_name, "tsv", sep="."), quote=F, row.names=F, col.names=T, sep="\t")


