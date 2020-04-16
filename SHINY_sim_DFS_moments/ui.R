library(shiny)
library(shinyWidgets)

fluidPage(
    titlePanel("Simulated DFS"),
    fluidRow(
        column(2,
               wellPanel(
                   h2("Add Migration"),

                   h3("Period 2"),
                   checkboxInput("mig_p2", NULL),
                   conditionalPanel(condition = "input.mig_p2 == true",
                   shinyWidgets::sliderTextInput("p2_m12_3", label = "pop12 <- pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p2_m3_12", label = "pop12 -> pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F)
                   ),
                   
                   h3("Period 3.1"),
                   checkboxInput("mig_p31", NULL),
                   conditionalPanel(condition = "input.mig_p31 == true",
                   shinyWidgets::sliderTextInput("p31_m23", label = "pop2 <- pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p31_m32", label = "pop2 -> pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p31_m13", label = "pop1 <- pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p31_m31", label = "pop1 -> pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p31_m12", label = "pop1 <- pop2", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p31_m21", label = "pop1 -> pop2", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F)
                   ),
                   
                   h3("Period 3.2"),
                   checkboxInput("mig_p32", NULL),
                   conditionalPanel(condition = "input.mig_p32 == true",
                   shinyWidgets::sliderTextInput("p32_m23", label = "pop2 <- pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p32_m32", label = "pop2 -> pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p32_m13", label = "pop1 <- pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p32_m31", label = "pop1 -> pop3", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p32_m12", label = "pop1 <- pop2", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F),
                   shinyWidgets::sliderTextInput("p32_m21", label = "pop1 -> pop2", choices = c(seq(0,1,0.1),seq(2,10,1)), selected=0, grid=F)
                   )
               )
        ),
        column(2,
               wellPanel(
                   h2("Modify Population Sizes (relative)"),
                   h3("Period 2"),
                   checkboxInput("mod_N_p2", NULL),
                   conditionalPanel(condition = "input.mod_N_p2 == true",
                   shinyWidgets::sliderTextInput("p2_N12", label = "Pop 1,2 ancestor", choices = c(0.1,0.5,1,5,10), selected=1, grid=F),
                   shinyWidgets::sliderTextInput("p2_N3", label = "Pop 3", choices = c(0.1,0.5,1,5,10), selected=1, grid=F)
                   ),
                   
                   h3("Period 3.1"),
                   checkboxInput("mod_N_p31", NULL),
                   conditionalPanel(condition = "input.mod_N_p31 == true",
                   shinyWidgets::sliderTextInput("p31_N1", label = "Pop 1", choices = c(0.1,0.5,1,5,10), selected=1, grid=F),
                   shinyWidgets::sliderTextInput("p31_N2", label = "Pop 2", choices = c(0.1,0.5,1,5,10), selected=1, grid=F),
                   shinyWidgets::sliderTextInput("p31_N3", label = "Pop 3", choices = c(0.1,0.5,1,5,10), selected=1, grid=F)
                   ),
                   
                   h3("Period 3.2"),
                   checkboxInput("mod_N_p32", NULL),
                   conditionalPanel(condition = "input.mod_N_p32 == true",
                   shinyWidgets::sliderTextInput("p32_N1", label = "Pop 1", choices = c(0.1,0.5,1,5,10), selected=1, grid=F),
                   shinyWidgets::sliderTextInput("p32_N2", label = "Pop 2", choices = c(0.1,0.5,1,5,10), selected=1, grid=F),
                   shinyWidgets::sliderTextInput("p32_N3", label = "Pop 3", choices = c(0.1,0.5,1,5,10), selected=1, grid=F)
                   )
               )
        ),
        column(2,
                     
               wellPanel(
                   h2("Times (2N generations)"),
                   shinyWidgets::sliderTextInput("p2_T", label = "Period 2", choices = c(0,seq(0.01,0.1,0.01),seq(0.2,1,0.1),seq(2,10,1)), selected=0.1, grid=F),
                   shinyWidgets::sliderTextInput("p31_T", label = "Period 3.1", choices = c(0,seq(0.01,0.1,0.01),seq(0.2,1,0.1),seq(2,10,1)), selected=0.1, grid=F),
                   shinyWidgets::sliderTextInput("p32_T", label = "Period 3.2", choices = c(0,seq(0.01,0.1,0.01),seq(0.2,1,0.1),seq(2,10,1)), selected=0.1, grid=F)
                ),
                
                wellPanel(
                   h2("Sampling"),
                   numericInput("n1_2", label = "Samples from pops 1 & 2", value = 30),
                   numericInput("n3", label = "Samples from pop 3", value = 10),
                   checkboxInput("swap13", label = "Swap pops 1 and 3")
                )
        ),
        column(6,
               plotOutput("model_plot"),
               actionButton("go", "Run diffusion approximation and plot DFS!"),
               plotOutput("DFS_plot"),
               checkboxInput("dcfs", "Plot dcfs"),
               conditionalPanel(condition = "input.dcfs == true",
               plotOutput("dcfs_plot")
               ),
               wellPanel(
                   textInput("prefix", label = "File prefix for saving", value = "DFS_simulation"),
                   actionButton("save_files", "Save plots and data")
                   )
        )
    )
)

 
