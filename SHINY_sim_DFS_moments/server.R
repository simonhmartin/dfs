source("DFS.R")

if.null <- function(x,y){
    if (is.null(x)==TRUE) return(y)
    return(x)
    }

make_3pop_sim_command <- function(params){
    paste("python3", "sim_SFS_moments.py",
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

plot.model <- function(params=NULL,
                       p1_T=1, p2_T=1, p31_T=1, p32_T=1,
                       p1_N=1,
                       p2_N12=1, p2_N3=1,
                       p31_N1=1, p31_N2=1, p31_N3=1,
                       p32_N1=1, p32_N2=1, p32_N3=1,
                       p2_m12_3=0, p2_m3_12=0,
                       p31_m12=0, p31_m21=0, p31_m13=0, p31_m31=0, p31_m23=0, p31_m32=0,
                       p32_m12=0, p32_m21=0, p32_m13=0, p32_m31=0, p32_m23=0, p32_m32=0,
                       label_periods=FALSE, arrow_col="black", reverse_populations=FALSE){
    
    try({p2_T=as.numeric(params$p2_T)})
    try({p31_T=as.numeric(params$p31_T)})
    try({p32_T=as.numeric(params$p32_T)})
    try({p2_N12=as.numeric(params$p2_N12)})
    try({p2_N3=as.numeric(params$p2_N3)})
    try({p31_N1=as.numeric(params$p31_N1)})
    try({p31_N2=as.numeric(params$p31_N2)})
    try({p31_N3=as.numeric(params$p31_N3)})
    try({p32_N1=as.numeric(params$p32_N1)})
    try({p32_N2=as.numeric(params$p32_N2)})
    try({p32_N3=as.numeric(params$p32_N3)})
    try({p2_m12_3=as.numeric(params$p2_m12_3)})
    try({p2_m3_12=as.numeric(params$p2_m3_12)})
    try({p31_m12=as.numeric(params$p31_m12)})
    try({p31_m21=as.numeric(params$p31_m21)})
    try({p31_m13=as.numeric(params$p31_m13)})
    try({p31_m31=as.numeric(params$p31_m31)})
    try({p31_m23=as.numeric(params$p31_m23)})
    try({p31_m32=as.numeric(params$p31_m32)})
    try({p32_m12=as.numeric(params$p32_m12)})
    try({p32_m21=as.numeric(params$p32_m21)})
    try({p32_m13=as.numeric(params$p32_m13)})
    try({p32_m31=as.numeric(params$p32_m31)})
    try({p32_m23=as.numeric(params$p32_m23)})
    try({p32_m32=as.numeric(params$p32_m32)})
    
    maxN <- max(p1_N, p2_N12, p2_N3, p31_N1, p31_N2, p31_N3, p32_N1, p32_N2, p32_N3)

    dist <- maxN+1

    p3_mid <- c(0, dist, 2*dist)

    p2_mid <- c(sum(p3_mid[1:2])/2, p3_mid[3])

    p1_mid <- sum(p2_mid)/2

    if (reverse_populations==FALSE){
        left <- min(c(p1_mid, p2_mid, p3_mid)) - 0.5*maxN
        right <- max(c(p1_mid, p2_mid, p3_mid)) + 0.5*maxN
        }
    else{
        right <- min(c(p1_mid, p2_mid, p3_mid)) - 0.5*maxN
        left <- max(c(p1_mid, p2_mid, p3_mid)) + 0.5*maxN
        }
    
    levels <- cumsum(rev(c(p1_T, p2_T, p31_T, p32_T)))
    top <- max(levels)
    
    plot(0, cex=0, ylim = c(0, top), xlim = c(left,right), bty="n", xaxt="n", ylab = "Generations ago (units of 2N)", xlab = "")

    rect(p3_mid - c(p32_N1, p32_N2, p32_N3)/2, 0, p3_mid + c(p32_N1, p32_N2, p32_N3)/2, levels[1], col="gray50", border=NA)
    rect(p3_mid - c(p31_N1, p31_N2, p31_N3)/2, levels[1], p3_mid + c(p31_N1, p31_N2, p31_N3)/2, levels[2], col="gray60", border=NA)
    rect(p2_mid - c(p2_N12, p2_N3)/2, levels[2], p2_mid + c(p2_N12, p2_N3)/2, levels[3], col="gray70", border=NA)
    rect(p1_mid - p1_N/2, levels[3], p1_mid + p1_N/2, levels[4], col="gray80", border=NA)
    
    segments(p3_mid[1] - p31_N1/2, levels[2], p3_mid[2] + p31_N2/2, levels[2], col="gray70")
    segments(p2_mid[1] - p2_N12/2, levels[3], p2_mid[2] + p2_N3/2, levels[3], col="gray80")
    
    if (p2_m12_3 > 0) arrows(p2_mid[2]-p2_N3/2, levels[2]+p2_T*.33, p2_mid[1]+p2_N12/2, levels[2]+p2_T*.33, lwd=p2_m12_3, col=arrow_col)
    if (p2_m3_12 > 0) arrows(p2_mid[1]+p2_N12/2, levels[2]+p2_T*.66, p2_mid[2]-p2_N3/2, levels[2]+p2_T*.66, lwd=p2_m3_12, col=arrow_col)
    
    if (p31_m23 > 0) arrows(p3_mid[3]-p31_N3/2, levels[1]+p31_T*.33, p3_mid[2]+p31_N2/2, levels[1]+p31_T*.33, lwd=p31_m23, col=arrow_col)
    if (p31_m32 > 0) arrows(p3_mid[2]+p31_N2/2, levels[1]+p31_T*.66, p3_mid[3]-p31_N3/2, levels[1]+p31_T*.66, lwd=p31_m32, col=arrow_col)
    
    if (p31_m13 > 0) arrows(p3_mid[3]-p31_N3/2, levels[1]+p31_T*.33, p3_mid[1]+p31_N1/2, levels[1]+p31_T*.33, lwd=p31_m13, col=arrow_col)
    if (p31_m31 > 0) arrows(p3_mid[1]+p31_N1/2, levels[1]+p31_T*.66, p3_mid[3]-p31_N3/2, levels[1]+p31_T*.66, lwd=p31_m31, col=arrow_col)

    if (p31_m12 > 0) arrows(p3_mid[2]-p31_N2/2, levels[1]+p31_T*.45, p3_mid[1]+p31_N1/2, levels[1]+p31_T*.45, lwd=p31_m12, col=arrow_col)
    if (p31_m21 > 0) arrows(p3_mid[1]+p31_N1/2, levels[1]+p31_T*.55, p3_mid[2]-p31_N2/2, levels[1]+p31_T*.55, lwd=p31_m21, col=arrow_col)
    
    if (p32_m23 > 0) arrows(p3_mid[3]-p32_N3/2, p32_T*.33, p3_mid[2]+p32_N2/2, p32_T*.33, lwd=p32_m23, col=arrow_col)
    if (p32_m32 > 0) arrows(p3_mid[2]+p32_N2/2, p32_T*.66, p3_mid[3]-p32_N3/2, p32_T*.66, lwd=p32_m32, col=arrow_col)
    
    if (p32_m13 > 0) arrows(p3_mid[3]-p32_N3/2, p32_T*.33, p3_mid[1]+p32_N1/2, p32_T*.33, lwd=p32_m13, col=arrow_col)
    if (p32_m31 > 0) arrows(p3_mid[1]+p32_N1/2, p32_T*.66, p3_mid[3]-p32_N3/2, p32_T*.66, lwd=p32_m31, col=arrow_col)
    
    if (p32_m12 > 0) arrows(p3_mid[2]-p32_N2/2, p32_T*.45, p3_mid[1]+p32_N1/2, p32_T*.45, lwd=p32_m12, col=arrow_col)
    if (p32_m21 > 0) arrows(p3_mid[1]+p32_N1/2, p32_T*.55, p3_mid[2]-p32_N2/2, p32_T*.55, lwd=p32_m21, col=arrow_col)
    
    if (reverse_populations == FALSE) mtext(1, at=p3_mid, text=1:3)
    else mtext(1, at=p3_mid, text=3:1)
    
    if(label_periods==TRUE){
        mtext(4,at=c(p32_T*.5, levels[1]+p31_T*.5, levels[2]+p2_T*.5),
            text=c("Period 3.2", "Period 3.1", "Period 2"), las=2)
        }
    }


run.simulation <- function(input){
    #convert population numbers for dadi. 3 becomes 1 and 1 becomes 3
    moments_params <- list(n1 = ifelse(input$swap13==FALSE, input$n3, input$n1_2),
                    n2 = input$n1_2,
                    n3 = ifelse(input$swap13==FALSE, input$n1_2, input$n3), #number of samples in each pop
                    
                    p2_T = input$p2_T, #length of two pop periods (in 2N generations)
                    p2_m21 = input$p2_m12_3, #first three pop period m23
                    p2_m12 = input$p2_m3_12, #first three pop period m32
                    p2_N1 = input$p2_N3, #first twp pop period N1
                    p2_N2 = input$p2_N12, #first twp pop period N1
                    
                    p31_T = input$p31_T, #first three pop period
                    p31_N1 = input$p31_N3, #first three pop period N1
                    p31_N2 = input$p31_N2, #first three pop period N2
                    p31_N3 = input$p31_N1, #first three pop period N2
                    p31_m21 = input$p31_m23, #first three pop period m23
                    p31_m12 = input$p31_m32, #first three pop period m32
                    p31_m31 = input$p31_m13, #first three pop period m31
                    p31_m13 = input$p31_m31, #first three pop period m13
                    p31_m23 = input$p31_m21, #first three pop period m23
                    p31_m32 = input$p31_m12, #first three pop period m32
                    
                    p32_T = input$p32_T, #second three pop period
                    p32_N1 = input$p32_N3, #second three pop period N1
                    p32_N2 = input$p32_N2, #second three pop period N2
                    p32_N3 = input$p32_N1, #second three pop period N2
                    p32_m21 = input$p32_m23, #second three pop period m23
                    p32_m12 = input$p32_m32, #second three pop period m32
                    p32_m31 = input$p32_m13, #second three pop period m31
                    p32_m13 = input$p32_m31, #second three pop period m13
                    p32_m23 = input$p32_m21, #second three pop period m23
                    p32_m32 = input$p32_m12 #second three pop period m32
    )
    
    command <- make_3pop_sim_command(moments_params)
    
    print(command)
        
    read.table(text=system(command, intern=TRUE,ignore.stderr=FALSE))
    }


shinyServer(
function(input, output, session){
    
    reactive_data <- reactiveValues()
    
    output$model_plot <- renderPlot({
        
        par(mar=c(1,15,1,15))
        
        plot.model(params=input,
                   p1_T=mean(as.numeric(input$p2_T), as.numeric(input$p31_T), as.numeric(input$p32_T))/2,
                   label_periods=TRUE, arrow_col="red", reverse_populations=input$swap13)
        
        }, height=300)
    
    
    output$DFS_plot <- renderPlot({
        
        get_SFS <- eventReactive(input$go, {run.simulation(input)})
        
        reactive_data$SFS <- get_SFS()
        
        par(mar=c(1.2,6,2,4))
        
        if (is.null(reactive_data$SFS) == FALSE){
            if (input$swap13 == FALSE){
                reactive_data$P1 <- 3
                reactive_data$P3 <- 1
                }
            else {
                reactive_data$P1 <- 1
                reactive_data$P3 <- 3
                }
            reactive_data$dfs_data <- get.DFS(base_counts=reactive_data$SFS[,c(reactive_data$P1,2,reactive_data$P3)],
                                                site_counts=reactive_data$SFS[,4],
                                                Ns=as.numeric(c(input$n1_2,input$n1_2,input$n3)))
            
            plotDFS(reactive_data$dfs_data$DFS, reactive_data$dfs_data$weights, col_D="red")
            
            }
        
        }, height=300)
    
    
    output$dcfs_plot <- renderPlot({
        if (is.null(reactive_data$SFS) == FALSE){
            if (input$dcfs == TRUE){
                reactive_data$dcfs <- get.dcfs(base_counts=reactive_data$SFS[,c(reactive_data$P1,2,reactive_data$P3)],
                                            site_counts=reactive_data$SFS[,4],
                                            Ns=as.numeric(c(input$n1_2,input$n1_2,input$n3)))
                par(mar=c(1,6,2,4))
                plot.dcfs(reactive_data$dcfs)
                }
            }
        }, height=300)

    
    observeEvent(input$save_files, {
        pdf(paste0(input$prefix,".model.pdf"), width=5,height=5)
        par(mar=c(1,5,1,1))
        plot.model(params=input, p1_T=mean(as.numeric(input$p2_T), as.numeric(input$p31_T), as.numeric(input$p32_T))/2, arrow_col="red", reverse_populations=input$swap13)
        dev.off()
    
        png(paste0(input$prefix,".model.png"), width = 1000, height = 800, res=300)
        par(mar=c(1,5,1,1))
        plot.model(params=input, p1_T=mean(as.numeric(input$p2_T), as.numeric(input$p31_T), as.numeric(input$p32_T))/2, arrow_col="red", reverse_populations=input$swap13)
        dev.off()
        
        pdf(paste0(input$prefix,".pdf"), width=10,height=4)
        par(mar=c(1.2,4,1,4))
        plotDFS(reactive_data$dfs_data$DFS, reactive_data$dfs_data$weights, col_D="red")
        dev.off()
        
        png(paste0(input$prefix,".png"), width = 2000, height = 800, res=300)
        par(mar=c(1.2,4,1,4))
        plotDFS(reactive_data$dfs_data$DFS, reactive_data$dfs_data$weights, col_D="red")
        dev.off()
        
        write.table(data.frame(D=round(reactive_data$dfs_data$DFS,4),
                               weight=round(reactive_data$dfs_data$weights,4)),
                    file=paste0(input$prefix,".csv"), sep=",",row.names=FALSE,quote=FALSE)
        
        write.table(reactiveValuesToList(input), file=paste0(input$prefix,".params.csv"), quote=F, row.names=F, sep=",")
        })

    }
    )

