library(shiny)
library(shinyalert)
library(shinyjs)
library(ProteoBayes)
library(reactable)
library(tidyverse)
library(grid)
`%notin%` <- Negate(`%in%`)

options(shiny.maxRequestSize = 50 * 1024^2) 

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "lumen", bg = "white", fg = "green"),

  # Application title
  tabsetPanel(
    tabPanel('Home',
             titlePanel('Welcome to ProteoBayes!'),
             
             useShinyjs(),
             
             fluidPage(column(width = 12), 
                      p("ProteoBayes is a Bayesian statistical framework for
                        differential proteomics analysis. It relies on mass 
                        spectrometry-based intensity measurements to be 
                        compared between biological conditions.",
                        align="justify", 
                        style = "color:black"),
                      p("More details can be found in the preprint by", 
                        tags$a("Chion & Leroy (2023).", href="https://arxiv.org/abs/2307.08975"),
                        align="justify", 
                        style = "color:black"),
                      p("The present Web-app is a code-free implementation of the
                        ProteoBayes R package, available on the",
                        tags$a("CRAN.", href="https://cran.r-project.org/web/packages/ProteoBayes/index.html"),
                        "Its development version can be found on",
                        tags$a(href="https://github.com/mariechion/ProteoBayes", "GitHub."),
                        align="justify", 
                        style = "color:black"),
                      p("A",
                        tags$a(href="https://mariechion.github.io/ProteoBayes/",
                               "documentation website"),
                        "is also provided to provides help and examples of its usage.",
                        align="justify",
                        style = "color:black"),
                      h3("Illustration of ProteoBayes with synthetic data",
                         align="justify"),
                      p("This tab allows you to try ProteoBayes with synthetic (simulated) 
                        datasets.",
                        align="justify",
                        style = "color:black"),
                      h3("ProteoBayes on a csv file",
                         align="justify"),
                      p("This tab allows you to use ProteoBayes with your own
                        data. Note that the submitted dataset needs to be a .csv 
                        file in long format with the following variables (column names):
                        'Peptide', 'Group', 'Sample' & 'Output' (for intensity values). 
                        See the previous tab with synthetic data for examples.",
                        align="justify",
                        style = "color:black"),
                      h3("Export results as table", 
                         align="justify"),
                      p("This tab allows you to export the results as a table.
                        Those results contain the comparison across all groups 
                        of their posterior means and 95% credible intervals 
                        computed for each peptide. A additional 'Distinct' column
                        indicates whether the groups are considered different or not
                        (this criterion is verified if credible intervals do not
                        overlap)." ,
                        align="justify",
                        style = "color:black"),
                      br(),
                      p("Any question? any feedback? Please send them to 
                        Marie Chion (mc2411{at}cam.ac.uk).")
             )
             
             
             ),
    
    
    tabPanel('Illustration of ProteoBayes with synthetic data',
             titlePanel('Compute mean differences from a synthetic dataset'),

             useShinyjs(),

             sidebarLayout(
               sidebarPanel(div(style = "overflow-y: auto; max-height: 400px;",  # Ajout des styles pour la barre de défilement
                                radioButtons("customizeSimuDb", "Modifying the dataset?", c('No','Yes')),
                                conditionalPanel(
                                  condition = "input.customizeSimuDb == 'Yes'",
                                  sliderInput("nb_peptide","Number of peptides", val = 5, min = 1, max = 1000, step = 1),
                                  sliderInput("nb_group","Number of groups/conditions", val = 5, min = 1, max = 20, step = 1),
                                  sliderInput("nb_sample","Number of samples per group", val = 5, min = 1, max = 100, step = 1),
                                  sliderInput("diff_group","Mean difference between consecutive groups", val = 3, min = 0, max = 10, step = 0.1)
                                ),
                                radioButtons("custom_hp_simu", "Adapting prior parameters?", c('No','Yes')),
                                conditionalPanel(
                                  condition = "input.custom_hp_simu == 'Yes'",
                                  sliderInput("mu_0", HTML("Prior value of &mu; coefficient"),
                                              value = 0, min = 0, max = 100, step = 0.1),
                                  sliderInput("lambda_0", HTML("Prior value of &lambda; coefficient"),
                                              value = 1, min = 0.1, max = 10, step = 0.1),
                                  sliderInput("alpha_0", HTML("Prior value of &alpha; coefficient"),
                                              value = 1, min = 0.1, max = 10, step = 0.1),
                                  sliderInput("beta_0", HTML("Prior value of &beta; coefficient"),
                                              value = 1, min = 0.1, max = 10, step = 0.1)
                                ),
                                uiOutput('simu_peptide_graph'),
                                uiOutput('simu_group1_graph'),
                                uiOutput('simu_group2_graph'),
                                actionButton("simu_and_plot_data", "Generate data and display results", icon=icon("chart-line")))),

               mainPanel(reactableOutput("tibble_simu"),
                         plotOutput("plot_simu")
               )
             )

    ),
    
    tabPanel('ProteoBayes on a csv file',
             titlePanel('Compute mean differences from a csv file'),

             useShinyjs(),

             sidebarLayout(
               sidebarPanel(div(style = "overflow-y: auto; max-height: 400px; width: 100%",  # Ajout des styles pour la barre de défilement
                                fileInput("upload",
                                          "Upload a csv file",
                                          accept = c('csv', 'comma-separated-values', ".csv"),
                                          buttonLabel = "Upload..."),
                                actionButton("display_data_csv", "Display data from csv", icon=icon("file")),
                                radioButtons("custom_hp_csv", "Adapting prior parameters?", c('No','Yes')),
                                conditionalPanel(
                                  condition = "input.custom_hp_csv == 'Yes'",
                                sliderInput("csv_mu_0", HTML("Prior value of &mu; coefficient"),
                                            value = 0, min = 0, max = 100, step = 0.1),
                                sliderInput("csv_lambda_0", HTML("Prior value of &lambda; coefficient"),
                                            value = 1, min = 0.1, max = 10, step = 0.1),
                                sliderInput("csv_alpha_0", HTML("Prior value of &alpha; coefficient"),
                                            value = 1, min = 0.1, max = 10, step = 0.1),
                                sliderInput("csv_beta_0", HTML("Prior value of &beta; coefficient"),
                                            value = 1, min = 0.1, max = 10, step = 0.1),
                                ),
                                uiOutput('csv_peptide_graph'),
                                uiOutput('csv_group1_graph'),
                                uiOutput('csv_group2_graph'),
                                actionButton("view_graph", "Display graph", icon=icon("chart-line")),
                                downloadButton("download_graphic", "Download Graphic", icon = icon("download")),
               )),

               mainPanel(uiOutput("error_panel1"),
                         uiOutput("warning_csv"),
                         uiOutput("error_colnames"),
                         uiOutput("error_format_colnames"),
                         reactableOutput("data_csv"),
                         plotOutput("plot"))
             )

    ),

    tabPanel('Export results as table',
             titlePanel('Export ProteoBayes results computed on your csv file'),
             div(style = "position: absolute; top: 100px; left: 20px;",
                 uiOutput("error_panel_export")),
             div(style = "position: absolute; top: 120px; left: 50px; width: 70%;",
                 reactableOutput("tibble_ProteoBayes")
             ),

             div(style = "position: absolute; top: 100px; right: 150px;",
                 downloadButton("download_data", "Download data",
                                icon = icon("download")))
    )
  )
)





# Define server logic
server <- function(input, output) {
  thematic::thematic_shiny()
  
  multi_imp <- reactive({
    if(input$multi_imp == 'Yes'){bool = TRUE} else {bool=FALSE}
    return(bool)
  })
  
  simulated_db <- reactive({
    ProteoBayes::simu_db(
      nb_peptide = input$nb_peptide,
      nb_group = input$nb_group,
      nb_sample = input$nb_sample,
      diff_group = input$diff_group,
    ) %>% 
      return()
  })
  
  observeEvent(input$simu_and_plot_data, {
  
    output$simu_peptide_graph = renderUI({
      selectInput('simu_peptide_graph', 'Peptide to analyse', 
                  choices = simulated_db()$Peptide %>% unique(), 
                  selected = unique(simulated_db()$Peptide)[1])
    })
    
    output$simu_group1_graph = renderUI({
      selectInput('simu_group1_graph', 'Group1 to compare', 
                  choices = unique(simulated_db()$Group),
                  selected = unique(simulated_db()$Group)[1])
    })
    
    output$simu_group2_graph = renderUI({
      selectInput('simu_group2_graph', 'Group2 to compare', 
                  choices = unique(simulated_db()$Group),
                  selected = unique(simulated_db()$Group)[2])
    })
    
    output$tibble_simu <- renderReactable(
      reactable(simulated_db(), defaultPageSize = 5))
    
    res_simu <- reactive({
      ProteoBayes::posterior_mean(
        data = simulated_db(),
        mu_0 = input$mu_0,
        lambda_0 = input$lambda_0,
        alpha_0 = input$alpha_0,
        beta_0 = input$beta_0) 
    })
    
    output$plot_simu <- renderPlot({
      set.seed(1)
       res_simu() %>% 
        ProteoBayes::sample_distrib(nb_sample = 1000) %>% 
        ProteoBayes::plot_distrib(
          group1 = input$simu_group1_graph,
          group2 = input$simu_group2_graph,
          peptide = input$simu_peptide_graph)
    })

  })  

  ########### Panel 2 ##############
  
  #data <- NULL
  data <- reactive({
    req(input$upload)

    ext <- tools::file_ext(input$upload$name)
    switch(
      ext,
      csv = vroom::vroom(input$upload$datapath, delim = NULL),
      validate("Invalid file; Please upload a .csv file")
    )
  })

  observeEvent(input$display_data_csv, {
    output$data_csv <- renderReactable(reactable(data(),defaultPageSize = 5))
  })
  
  output$csv_peptide_graph = renderUI({
    selectInput('csv_peptide_graph', 'Peptide to analyse', 
                choices = data()$Peptide %>% unique(), 
                selected = unique(data()$Peptide)[1])
  })
  
  output$csv_group1_graph = renderUI({
    selectInput('csv_group1_graph', 'Group1 to compare', 
                choices = unique(data()$Group),
                selected = unique(data()$Group)[1])
  })
  
  output$csv_group2_graph = renderUI({
    selectInput('csv_group2_graph', 'Group2 to compare', 
                choices = unique(data()$Group),
                selected = unique(data()$Group)[2])
  })
  
  res <- reactive({
    ProteoBayes::posterior_mean(
      data = data(),
      mu_0 = input$csv_mu_0,
      lambda_0 = input$csv_lambda_0,
      alpha_0 = input$csv_alpha_0,
      beta_0 = input$csv_beta_0)
  })

  reac_button <- eventReactive(input$view_graph, {
    set.seed(1)
    res() %>% 
      ProteoBayes::sample_distrib() %>%  
      ProteoBayes::plot_distrib(
        group1 = input$csv_group1_graph,
        group2 = input$csv_group2_graph,
        peptide = input$csv_peptide_graph) %>% 
      return()
  })
  
  output$plot <- renderPlot({
    reac_button()
  })

  output$download_graphic <- downloadHandler(
    filename = "ProteoBayes_plot.png",
    content = function(file) {
      plot_csv <- reac_button()
      if (!is.null(plot_csv)) {
        ggsave(
          file,
          plot = plot_csv,
          device = "png",
          width = unit(800 / 96, "in"),
          height = unit(600 / 96, "in")
        )
      }
    }
  )


  ############# Panel 3 ###############

  output$tibble_ProteoBayes <- renderReactable(
    reactable(res() %>% identify_diff()))

  output$download_data <- downloadHandler(
    filename = function() {
      paste("ProteoBayes_data.csv")
    },
    content = function(file) {
      vroom::vroom_write(res(), file, delim = ',')
    }
  )

}

# Run the application
options(shiny.error = function() {
  return(NULL)
})

shinyApp(ui = ui, server = server)


