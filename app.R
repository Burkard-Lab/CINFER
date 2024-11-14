#This is application is meant to provide others the means to infer rates of mis-segregation in their own samples using the analytical process outlined in Lynch et al. 2022 (DOI: 10.7554/eLife.69799).

# Packages ----------------------------------------------------------------
library(shiny)
library(htmltools)
library(bslib)
library(bsicons)
library(thematic)
library(shinyWidgets)
library(shinycssloaders)
library(plyr)
library(tidyverse)
library(ggpubr)
library(ggExtra)
library(abc)
library(phyloTop)
library(dbplyr)
library(pheatmap)
library(DT)
library(reshape2)
library(uwot)
library(nlrx)
library(progressr)
library(splitstackshape)
library(gtools)
library(shinyjs)
library(grid)
library(odbc)
library(DBI)
library(pool)

# Functions ---------------------------------------------------------------

linebreaks <- function(n){HTML(strrep(br(), n))}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

source("src/CINFERPostPred.R")

theme_cinfer <- function(){
  font <- "Arial"
  coord_cartesian(clip = "off")
  theme_classic() %+replace%
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 1, color = "grey50", fill = "transparent"),
      axis.line.y = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks = element_line(linewidth = 0.5, color = "grey50"),
      axis.line = element_blank(),
      axis.text = element_text(size = 10, color = "grey50"),
      axis.title = element_text(size = 10, face = "bold", color = "grey50"),
      strip.background = element_blank(),
      strip.text = element_text(size = 10)
    )
}

chrs <- paste("chr", c(seq(1,22), "X", "Y"), sep = "")
chrarms <- paste("chr", rep(c(seq(1,22), "X", "Y"), each = 2), c("p","q"), sep = "")

headerCallbackRemoveHeaderFooter <- c(
  "function(thead, data, start, end, display){",
  "  $('th', thead).css('display', 'none');",
  "}"
)

# Database ----------------------------------------------------------------
con <- dbConnect(RSQLite::SQLite(),
                dbname = "res/CINFERbasev1.0.sqlite")
results_db <- tbl(con, "Exp_Dip_Abun_sumstats_small")

# pool <- dbPool(
#   drv = RMariaDB::MariaDB(), 
#   user = Sys.getenv("AWS_USER"),
#   password = Sys.getenv("AWS_PASS"),
#   dbname = Sys.getenv("AWS_DBNAME"),
#   host = Sys.getenv("AWS_HOST")
# )

onStop(function() {
  dbDisconnect(con)
})

# User Interface ----------------------------------------------------------
ui <- fluidPage(
   useShinyjs(),
   tags$style("
              body {

    zoom: 0.8; /* Other non-webkit browsers */
    zoom: 80%; /* Webkit browsers */
}
              "),
  tags$head(
    HTML("<title>CINFER</title>"),
    tags$link(rel="shortcut icon", type="image/x-icon", href="res/favicon.ico")),
  theme = bslib::bs_theme(bootswatch = "simplex"),
  chooseSliderSkin("Flat"),
  div(titlePanel(title = div(strong(h1('CINFER', style="margin: 0; font-weight: bold; color: rgba(219, 46, 46, 0.8);")), h6('scDNA-based Estimation of Mis-segregation Rates', style="margin: 0;"))), style = "position:sticky;  position: -webkit-sticky; /* Safari */;top: 32px; min-height: 128px;"),
  sidebarLayout(

## Side Bar ----------------------------------------------------------------
      sidebarPanel(style = "overflow-y:scroll;position:fixed; max-height: 90%; width: 24%; max-width: 24%;",
                   strong("1 - Input Data"),
                   fileInput("file", "Upload your .csv file", buttonLabel = list(icon("glyphicon glyphicon-upload", lib="glyphicon"), "Upload"), accept = c(".csv")),
                   selectizeInput(
                     inputId = "testfile",
                     label = "... or select sample data", 
                     choices = list('', "Laks & McPherson et al, 2019" = setNames(c(
                       "ov2295",
                       "ov2295.clone1",
                        "ov2295.clone2",
                        "ov2295.clone3",
                        "ov2295.clone4"
                       ), 
                       c(
                         "OV2295 (all)",
                         "OV2295 (clone 1)",
                         "OV2295 (clone 2)",
                         "OV2295 (clone 3)",
                         "OV2295 (clone 4)"
                       )),
                      "Lynch et al, 2023" = setNames(
                        c("cal51","mcf10a"), 
                        c("CAL51","MCF10A"))),
                     multiple = F
                     ),
                   strong("2 - Select Priors"),
                   sliderInput("rateprior", h6("Mis-segregation Rate (MDD)"), min = 0, max = 4.6, value = c(0,2.3)),
                   sliderInput("selectionprior", h6("Selection Pressure (S)"), min = 0, max = 100, value = c(0,50)),
                   sliderInput("timeprior", h6("Time Scale (steps)"), min = 0, max = 100, value = c(0,50)),
                   #radioButtons("selmodel", "Selection Model", choices = c("Abundance", "Driver", "Hybrid", "Neutral"), selected = "Abundance"),
                   #radioButtons("growmodel", "Growth Model", choiceNames = c("Exponential Pseudo-Moran", "Constant Wright-Fisher"), choiceValues = c("exponential", "constant_wrightfisher"), selected = "exponential"),
                   pickerInput(
                     inputId = "selectedsumstats",
                     label = "Choose Summary Statistics", 
                     choices = list("Karyotype Diversity" = setNames(c("AvgAneuploidy", "NormMKV", "NormUniqueClones"), c("Aneuploidy", "MKV", "Unique Clones")),
                                    "Phylogenetic Topology" = setNames(c("SackinNorm", "CollessNorm", "CollessPermute", "NormCherries"), c("Sackin (norm)", "Colless (norm)", "Colless (permuted)", "Cherries"))),                     
                     selected = c("AvgAneuploidy", "NormMKV", "CollessPermute"),
                     multiple = TRUE,
                     options = list(
                       'selected-text-format' = "count > 3",
                       title = "Please choose at least 2"
                     )
                   ),
                   actionButton("parameterpull", "Pull Priors", icon=icon("glyphicon glyphicon-filter", lib="glyphicon")),
                   linebreaks(2),
                   strong("3 - Infer parameters"),
                   linebreaks(2),
                   radioButtons("abcmethod", "ABC Method", choiceNames = c("Rejection", "Rejection + local linear regression", "Rejection + ridge regression", " Rejection + neural net"), choiceValues = c("rejection", "loclinear", "ridge", "neuralnet")),
                   numericInput("tol", h6("Tolerance rate"), value = 0.05, min = 0.0001, max = 1),
                   #checkboxGroupInput("opts", h6("Other Options"), choices = c("Chromosome breakage? (P=0.5)")),
                   strong("4 - Run & Check ABC"),
                   linebreaks(1),
                   actionButton("abcrun", "Run", icon=icon("glyphicon glyphicon-equalizer", lib="glyphicon")),
                   actionButton("runchecks", "Check", icon=icon("glyphicon glyphicon-search", lib="glyphicon")),
                   linebreaks(2),
                   strong("5 - Report"),
                   prettyCheckboxGroup(
                     inputId = "reportitems",
                     label = "",
                     inline = T,
                     selected = "Stats",
                     choices = c("Stats", "Plots", "Data")
                   ),
                   downloadButton("report", "Download", icon=icon("glyphicon glyphicon-download", lib="glyphicon")),
        width = 3
      ),

## Main Panel --------------------------------------------------------------
      mainPanel(width = 9,
         fluidRow(
           layout_column_wrap(
              width = 1/2,
              style = css(grid_template_columns = "2fr 1fr"),
              card(
                 card_header(
                    class = "align-text-top fs-6",
                    strong("Input - Copy Number Profiles")
                 ),
                 card_body_fill(
                    plotOutput("cnheatmap", fill = T) %>% 
                      withSpinner(type = 6, color = "#e66e6e", size = 0.5)
                    ),
                 class="shadow-sm p-0 mb-3 bg-body rounded",
                 full_screen = T
                 ),
              card(
                 card_header(
                    class = "align-text-top fs-6",
                    strong("Input - Summary Statistics")
                 ),
                 card_body_fill(
                    dataTableOutput("inputstats") %>% 
                      withSpinner(type = 6, color = "#e66e6e", size = 0.5)
                 ),
                 class="shadow-sm p-0 mb-3 bg-body rounded"
              ),
              card(
                 card_header(
                    class = "align-text-top fs-6",
                    strong("Prior - Data")
                 ),
                 card_body_fill(
                    dataTableOutput("priorstats") %>% 
                      withSpinner(type = 6, color = "#e66e6e", size = 0.5)
                 ),
                 class="shadow-sm p-0 mb-3 bg-body rounded"
              ),
              card(
                 card_header(
                    class = "align-text-top fs-6",
                    strong("Prior - Summary Space")
                 ),
                 card_body_fill(
                    plotOutput("umap", fill = T) %>% 
                      withSpinner(type = 6, color = "#e66e6e", size = 0.5)
                 ),
                 class="shadow-sm p-0 mb-3 bg-body rounded"
              )
           )
           ),
         fluidRow(
            layout_column_wrap(
               width = 1/2,
               style = css(grid_template_columns = "1fr 2fr"),
               card(
                  card_header(
                     class = "align-text-top fs-6",
                     strong("Infer - Posterior Distribution")
                  ),
                  card_body_fill(
                  plotOutput("jointdists", fill = T) %>% withSpinner(type = 6, color = "#e66e6e", size = 0.5),
                  ),
                  class="shadow-sm p-0 mb-3 bg-body rounded"
               ),
               card(
                  card_header(
                     class = "align-text-top fs-6",
                     strong("Infer - Posterior Data")
                  ),
                  card_body_fill(
                    dataTableOutput("posteriorstats") %>% 
                      withSpinner(type = 6, color = "#e66e6e", size = 0.5)
                  ),
                  class="shadow-sm p-0 mb-3 bg-body rounded"
               )
            ),
            layout_column_wrap(
              width = 1/3,
              uiOutput("rateBox"),
              uiOutput("selBox"),
              uiOutput("timeBox")
            ),
            layout_column_wrap(
               width = 1/2,
               card(
                  card_header(
                     class = "align-text-top fs-6",
                     strong("Check - Posterior Prediction")
                  ),
                  card_body(
                    plotOutput("postpred") %>% 
                      withSpinner(type = 6, color = "#e66e6e", size = 0.5)
                  ),
                  class="shadow-sm p-0 mb-3 bg-body rounded"
               ),
               card(
                  card_header(
                     class = "align-text-top fs-6",
                     strong("Check - Tolerance Threshold")
                  ),
                  plotOutput("tolcv") %>% withSpinner(type = 6, color = "#e66e6e", size = 0.5),
                  class="shadow-sm p-0 mb-3 bg-body rounded"
               )
            )
         ),
        hr(),
        fluidRow(
           layout_column_wrap(
              width = 1,
              style = css(grid_template_columns = "4fr 1fr"),
              card(
                    card_header(strong("About"),
                                class = "align-text-top fs-6"),
              card_body(
                    p("CINFER is used to infer the rate of chromosome mis-segregation occuring in a sample directly from scDNAseq data, accounting for ongoing karyotype selection. CINFER uses an ", a(href = "https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf", "R implementation of Approximate Bayesian Computation"), " to make these inferences from a computational model of CIN developed in ", a(href="https://ccl.northwestern.edu/netlogo/", "NetLogo."), "We have shown that this method can recall experimentally observed mis-segregation rates as well as those most likely to have resulted in the karyotype diversity observed in patient-derived clinical samples.")
                        ),
              class="shadow-sm p-0 mb-3 bg-body rounded"
                 ),
              card(
                 card_header(strong("Bring me to..."),
                             class = "align-text-top fs-6"),
                 card_body_fill(
                 a(strong("the paper"), href="https://elifesciences.org/articles/69799"),
                 a(strong("the code"), href="https://github.com/Burkard-Lab/CINFER"),
                 a(strong("the Burkard Lab"), href="https://github.com/Burkard-Lab")
                 ),
                 class="shadow-sm p-0 mb-3 bg-body rounded"
              )
              # ,
              # card(
              #    card_header(strong("FAQ"),
              #                class = "align-text-top fs-6"),
              #    card_body(
              #      h6(strong("What is MDD?")),
              #      p("Mis-segregations per Diploid Division (MDD) is a standardized unit to for comparing chromosome mis-segregation rates across a number of different factors including ploidy, penetrance, and mis-segregations per abnormal division."),
              #       h6(strong("How should I choose my prior parameters?")),
              #       p(""),
              #       h6(strong("How do I interpret my results?")),
              #       h6(strong("What is ABC?")),
              #       p("Approximate Bayesian Computation is a method of statistical inference that determines the most likely set of biological parameters resulting in a given set of summary statistics based upon a simulated measurements of the system.")
              #       ),
              #    class="shadow-sm p-0 mb-3 bg-body rounded"
              #    )
           ),
              )
           )
      )
    )

# Server ------------------------------------------------------------------
server <- function(input, output) {
    thematic::thematic_shiny()
    shinyjs::disable("report")
  
  # Garbage Collection -----------------------------------------------------
    observeEvent(input$parameterpull, {
      req(input$testfile)
      gc()
    }) 
    
    observeEvent(input$abcrun, {
      req(input$testfile)
      gc()
    }) 
    
    observeEvent(input$runchecks, {
      req(input$testfile)
      gc()
    }) 

  # Input -------------------------------------------------------------------
    rv <- reactiveValues()
    
    observeEvent(input$file, {
      req(input$file)
      rv$inputfile = input$file$datapath
      reset("testfile")
    }) 
    
    observeEvent(input$testfile, {
      req(input$testfile)
      enable("testfile")
      rv$inputfile = paste("dat/",input$testfile,".csv",sep = "")
      reset("file")
    }) 
    
    inputcnv <- reactive({
      req(rv$inputfile)
      as.matrix(read.csv(rv$inputfile, header = T))
    })
    
    plot_cnheatmap <- reactive({
      pheatmap(inputcnv(),
                       cluster_cols = F,
                       border_color = "transparent",
                       treeheight_row = 100,
                       number_color = "grey50"
                       )
    })
    
    output$cnheatmap <- renderPlot({
      req(rv$inputfile)
      plot_cnheatmap()
    })
    
    inputstats <- reactive({
      req(rv$inputfile)
      chrMat <- inputcnv()
      clusMat <- hclust(dist(chrMat))
      
      colless.permute <- vector()
      for(i in 1:100){
        chrMat.permute <- apply(chrMat, 2, permute)
        chrClust.permute <- hclust(dist(chrMat.permute))
        colless.permute[i] <- colless.phylo(as.phylo(chrClust.permute), normalise = T)
      }
      modekaryo <- apply(chrMat, 2, getmode)
      chrMatsubclonal <- sweep(chrMat, 2, modekaryo, "-")
      data.frame(
                 AvgAneuploidy = mean(apply(chrMatsubclonal, 1, var, na.rm = T)), 
                 NormMKV = mean(apply(chrMat, 2, var)) / mean(apply(chrMat, 1, mean)), 
                 NormUniqueClones = nrow(unique(chrMat)) / nrow(chrMat),
                 SackinNorm = sackin.phylo(as.phylo(clusMat), normalise = TRUE), 
                 CollessNorm = colless.phylo(as.phylo(clusMat), normalise = TRUE),
                 CollessPermute = mean(colless.permute),
                 NormCherries = cherries(as.phylo(clusMat)) / nrow(chrMat) 
      )
    })
    
    output$inputstats <- renderDataTable(melt(inputstats()),
                                         options = list(
                                           dom = "t",
                                           paging = F,
                                           ordering = F,
                                           searching = F,
                                           headerCallback = JS(headerCallbackRemoveHeaderFooter)
                                         ),
                                         selection = 'none',
                                         callback = JS(
                                           "$('table.dataTable.no-footer').css('border-bottom', 'none');"
                                         ),
                                         class = 'row-border',
                                         escape = FALSE,
                                         rownames = FALSE,
                                         filter = "none")
    
    datsel <- eventReactive(input$parameterpull, {
      req(rv$inputfile, unlist(input$selectedsumstats) >= 2)
      
      #results_db <- tbl(pool, "Exp_Dip_Abun_sumstats_small")
      
      data.frame(results_db %>%
                   #filter(growth_model == !!input$growmodel) %>%
                   #filter(model == !!input$selmodel) %>%
                   filter(
                       between(Rate,!!input$rateprior[1], !!input$rateprior[2]) &
                       between(Pressure, !!input$selectionprior[1], !!input$selectionprior[2]) &
                       between(Step, !!input$timeprior[1], !!input$timeprior[2])) %>%
                   dplyr::select(c(Run, Step, Pressure, Rate, input$selectedsumstats)) %>%
                   collect())
    })
    
    output$priorstats <- renderDataTable(datsel(),
                                         options = list(
                                           pageLength = 6,
                                           dom = "tp"
                                         )
    )
    
    abctarget <- reactive({
      data.frame(inputstats()[,input$selectedsumstats])
    })
    
    abcstats <- reactive({
      data.frame(datsel()[,input$selectedsumstats])
    })
    
    abcpriors <- reactive({
      data.frame(datsel()[,c("Rate", "Pressure", "Step")])
    })
    
    abcresults <- eventReactive(input$abcrun,{
      req(rv$inputfile)
      res <- abc(target=abctarget(),
                 param=abcpriors(),
                 sumstat=abcstats(),
                 tol=input$tol,
                 method=input$abcmethod)
      return(res)
    })
    
    posteriordata <- eventReactive(input$abcrun, {
      abcdat <- data.frame(abcresults()[["unadj.values"]],abcstats()[which(abcresults()$region == TRUE),])
      abcdat$Rate <- abcdat$Rate * (46 / mean(apply(inputcnv(), 1, sum))) # Conversion to MDD
      shinyjs::enable("report")
      abcdat
    })
    
    output$rateBox <- renderUI({
      req(posteriordata())
      value <- round(median(data.frame(posteriordata())$Rate),3)
      stdev <- round(sd(data.frame(posteriordata())$Rate),3)
      card(
        bs_icon("speedometer2", size = "4em"),
        br(),
        h5("Mis-segregation Rate (MDD)"),
        br(),
        strong(h1(value)),
        br(),
        strong(h4(paste("±", stdev))),
        class = "shadow-sm p-0 mb-3 bg-dark rounded",
        style = "text-align: center; line-height: 0.5;"
      )
    })
    
    output$selBox <- renderUI({
      req(posteriordata())
      value <- round(median(data.frame(posteriordata())$Pressure),3)
      stdev <- round(sd(data.frame(posteriordata())$Pressure),3)
      card(
        bs_icon("stoplights", size = "4em"),
        br(),
        h5("Selection Pressure (S)"),
        br(),
        strong(h1(value)),
        br(),
        strong(h4(paste("±", stdev))),
        class = "shadow-sm p-0 mb-3 bg-dark rounded",
        style = "text-align: center; line-height: 0.5;"
      )
    })
    
    output$timeBox <- renderUI({
      req(posteriordata())
      value <- round(median(data.frame(posteriordata())$Step),3)
      stdev <- round(sd(data.frame(posteriordata())$Step),3)
      card(
        bs_icon("hourglass-split", size = "4em"),
        br(),
        h5("Time Steps"),
        br(),
        strong(h1(value)),
        br(),
        strong(h4(paste("±", stdev))),
        class = "shadow-sm p-0 mb-3 bg-dark rounded",
        style = "text-align: center; line-height: 0.5;"
      )
    })
    
    output$posteriorstats <- renderDataTable(posteriordata(),
                                             options = list(
                                               pageLength = 6,
                                               dom = "tp"
                                             )
                                             )
    
    plot_jointdists <- reactive({
      ggplot(data.frame(abcresults()[["unadj.values"]]), aes(Rate,Pressure)) + 
        stat_density_2d_filled(contour_var = "ndensity") + 
        geom_vline(aes(xintercept = median(Rate)), color = "black", linetype = 1) + 
        geom_hline(aes(yintercept = median(Pressure)), color = "black", linetype = 1) + 
        geom_vline(aes(xintercept = median(Rate)), color = "white", linetype = 2) + 
        geom_hline(aes(yintercept = median(Pressure)), color = "white", linetype = 2) + 
        theme_cinfer() + 
        theme(
           legend.position = "none",
           axis.title = element_text(size = 14)) +
        scale_x_continuous(expand = c(0,0)) + 
        scale_y_continuous(expand = c(0,0)) + 
          scale_fill_manual(values = colorRampPalette(c("transparent", "skyblue", "yellowgreen", "orange", "red", "darkred"))(10)) + 
          labs(x = "MDD", 
               y = "Pressure (S)")
    })
    
    output$jointdists <- renderPlot({
      req(rv$inputfile)
      plot_jointdists()
    })
    
    tolcvdata <- eventReactive(input$runchecks, {
       req(rv$inputfile)
       cvdata <- cv4abc(param = abcpriors()[,"Rate"], sumstat = abcstats(), nval = 3, tols = c(input$tol), method = input$abcmethod, statistic = "mean")
       cvdata.melt <- melt(data.frame("true" = cvdata$true[,1], "tol" = cvdata$estim[[1]]), id.vars = "true")
       return(cvdata.melt)
       })
    
    
    plot_tolcv <- reactive({
       ggplot(tolcvdata()) +
          geom_abline(slope = 1, color = "grey50") +
          geom_point(aes(x=true, y = value), size = 3, shape = 21, fill = "red") +
          theme_cinfer() +
          theme(
             legend.title = element_blank(),
             legend.background = element_blank(),
             legend.text = element_text(color = "grey50"),
             axis.title = element_text(size = 14)
          ) +
          #scale_fill_brewer(palette = "Reds") +
          labs(
             x = "True",
             y = "Estimate"
          )
    })

    output$tolcv <- renderPlot({
      req(rv$inputfile)
      plot_tolcv()
    })
    
     umap_data <- eventReactive(input$parameterpull, {
        req(rv$inputfile)
        unique_data <- unique.data.frame(abcstats())
        if(nrow(unique_data) > 2500){
           unique_data <- unique_data[sample(rownames(unique_data), 2500, replace = F),]
        }
        unique_data$group = "Sim"
        unique_data <- rbind(unique_data, data.frame(abctarget(), group = "Expt"))
        output_data <- data.frame(umap(unique_data, fast_sgd = T, n_neighbors = 5, bandwidth = 0.33), group = unique_data$group)
        output_data$group <- factor(output_data$group, levels = c("Sim", "Expt"))
        return(output_data)
     })
    
     plot_umap <- reactive({
       req(rv$inputfile)
       ggplot(umap_data(), aes(X1,X2)) +
         geom_point(aes(fill = group, size = group, color = group, shape = group)) +
         theme_cinfer() +
         theme(
           legend.position = "none",
           axis.title = element_text(size = 14)
         ) +
         scale_color_manual(values = c("grey50", "black")) +
         scale_fill_manual(values = c("grey50", "red")) + 
         scale_size_manual(values = c(2,4)) + 
         scale_shape_manual(values = c(21, 23)) + 
         labs(
           x = "UMAP1",
           y = "UMAP2"
         )
     })
     
     output$umap <- renderPlot({
        req(rv$inputfile)
        plot_umap()
       })
     
     postpred.data <- eventReactive(input$runchecks, {
        req(posteriordata())
        CINFERPostPred(#g = input$growmodel, 
                       #m = input$selmodel, 
                       r = median(data.frame(posteriordata())$Rate), 
                       s = median(data.frame(posteriordata())$Pressure), 
                       t = median(data.frame(posteriordata())$Step))
     })
     
     plot_postpred <- reactive({
       inputdata <- melt(data.frame(inputstats())[,unlist(input$selectedsumstats)])
       names(inputdata) <-  c("name", "value")
       
        ggplot(postpred.data() %>% filter(name %in% c(unlist(input$selectedsumstats))), aes(x = factor(s), y = value)) +
          facet_wrap(.~name, scales = "free", strip.position = "left") +
          geom_hline(data = inputdata, aes(yintercept = value*0.5), linetype = 2, color = "black") +
          geom_hline(data = inputdata, aes(yintercept = value*0.75), linetype = 2, color = "black") +
          geom_hline(data = inputdata, aes(yintercept = value), linetype = 1, color = "black") +
          geom_hline(data = inputdata, aes(yintercept = value*1.25), linetype = 2, color = "black") +
          geom_hline(data = inputdata, aes(yintercept = value*1.5), linetype = 2, color = "black") +
          stat_summary(fun = "mean", geom = "crossbar", width = 0.25) +
          geom_jitter(aes(fill = factor(s), color = factor(s)), size = 3, shape = 21, width  = 0.1) +
          theme_cinfer() + 
          theme(
              axis.title = element_text(size = 14),
              strip.placement = "outside",
              strip.text = element_text(size = 14, color = "grey50", face = "bold"),
              legend.position = "none"
          ) + 
          scale_fill_manual(values = c("grey50", "red")) +
          scale_color_manual(values = c("grey50", "black")) +
          labs(x = "Selection Pressure (S)",
               y = "")
     })
     
     output$postpred <- renderPlot({
       req(rv$inputfile, unlist(input$selectedsumstats) >= 2)
       plot_postpred()
     })
     
     
     output$report <- downloadHandler(
       filename = function() {
         paste(gsub(".csv", "", basename(file.path(rv$inputfile))), ".zip", sep = "")
       },
       content = function(file) {
         file_prefix = gsub(".csv", "", basename(file.path(rv$inputfile)))
         
         tmpdir = file.path(tempdir(), as.integer(Sys.time()))
         dir.create(tmpdir)
         
         if("Stats" %in% input$reportitems){
           write.csv(  
           data.frame(
               StatNames = paste(names(inputstats()), collapse = " "),
               SampleStats = paste(inputstats(), collapse = " "),
               SelectedStats = paste(input$selectedsumstats, collapse = " "),
                 Rate = paste(sort(unique(datsel()$Rate)), collapse = " "),
                 Pressure = paste(sort(unique(datsel()$Pressure)), collapse = " "),
                 Step = paste(sort(unique(datsel()$Step)), collapse = " "),
                 Rate_Min = min(posteriordata()$Rate),
                 Rate_Max = max(posteriordata()$Rate),
                 Rate_Mean = mean(posteriordata()$Rate),
                 Rate_Median = median(posteriordata()$Rate),
                 Rate_SD = sd(posteriordata()$Rate),
                 Rate_SE = plotrix::std.error(posteriordata()$Rate),
                 Pressure_Min = min(posteriordata()$Pressure),
                 Pressure_Max = max(posteriordata()$Pressure),
                 Pressure_Mean = mean(posteriordata()$Pressure),
                 Pressure_Median = median(posteriordata()$Pressure),
                 Pressure_SD = sd(posteriordata()$Pressure),
                 Pressure_SE = plotrix::std.error(posteriordata()$Pressure),
                 Step_Min = min(posteriordata()$Step),
                 Step_Max = max(posteriordata()$Step),
                 Step_Mean = mean(posteriordata()$Step),
                 Step_Median = median(posteriordata()$Step),
                 Step_SD = sd(posteriordata()$Step),
                 Step_SE = plotrix::std.error(posteriordata()$Step)
               ),
           file.path(tmpdir, paste(file_prefix, "Stats.csv", sep = "_")), row.names = F)
         }
         
         if("Plots" %in% input$reportitems){
           # Heatmap
           pdf(file = file.path(tmpdir, paste(file_prefix, "CopyHeatmap.pdf", sep = "_")))
           grid::grid.newpage()
           grid::grid.draw(plot_cnheatmap())
           dev.off()
           
           # Summary Space
           ggsave(plot = plot_umap(), file = file.path(tmpdir, paste(file_prefix, "SummarySpace.pdf", sep = "_")))
           
           # Posterior Distribution
           ggsave(plot = plot_jointdists(), file = file.path(tmpdir, paste(file_prefix, "PosteriorDistribution.pdf", sep = "_")))
           
           # Posterior Prediction and Tolerance Threshold (if generated)
           if(isTruthy(input$runchecks)){
             ggsave(plot = plot_postpred(), file = file.path(tmpdir, paste(file_prefix, "PosteriorPrediction.pdf", sep = "_")))
             ggsave(plot = plot_tolcv(), file = file.path(tmpdir, paste(file_prefix, "TolCV.pdf", sep = "_")))
           }
         }
         
         if("Data" %in% input$reportitems){
           write.csv(posteriordata(), file.path(tmpdir, paste(file_prefix, "PosteriorData.csv", sep = "_")), row.names = FALSE)
         }
         
         zip::zip(
           zipfile = file,
           files = dir(tmpdir),
           root = tmpdir
         )
       }
     )
}

# Run the application 
shinyApp(ui = ui, server = server)
