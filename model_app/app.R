library(shiny)
library(ggplot2)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Adaptive Sequencing Enrichment by Composition Model"),
  fluidRow(
    column(4,
  p("Adaptive sampling is a method of software-controlled enrichment unique to nanopore sequencing platforms recently implemented in Oxford Nanopore’s control software. 
      By examining the first few hundred bases of a DNA molecule as it passes through a pore, software can determine if the molecule is sufficiently interesting to sequence in its entirety. 
      If not, the molecule is ejected from the pore by reversing the voltage across it, freeing the pore for a new molecule. User supplied sequences define the targets to be sequenced or ejected.
      We developed a model to predict the enrichment-by-composition based on experimental parameters, and provided this web app so that researchers can explore the potential of adaptive sequencing for their own experiments. 
      Adjust the parameters below to see the effect it has on the model."),

      wellPanel(
        textInput(inputId = "read_length", label = "Average read Length:", value = "5000"),
        sliderInput(inputId = "species_abundance", label = "Species Abundance (%)", min = 0,  max = 100, value = 20, step = 0.1),
        sliderInput(inputId = "sequencing_speed", label = "Sequencing Speed (bps)", min = 0,  max = 1000, value = 420),
        sliderInput(inputId = "decision_time", label = "Decision Time (s)", min = 0,  max = 5, value = 1, step = 0.1),
        sliderInput(inputId = "capture_time", label = "Capture Time (s)", min = 0,  max = 5, value = 0.5, step = 0.1),
      ),
      p("Paper: S Martin, D Heavens, Y Lan, S Horsfield, MD Clark, RM Leggett.", strong("Nanopore adaptive sampling: a tool for enrichment of low abundance species in metagenomic sample."), em("bioRxiv 2021.05.07.443191"), "; doi: ", a("https://doi.org/10.1101/2021.05.07.443191", href="https://doi.org/10.1101/2021.05.07.443191"))
    ),

    column(8,
        plotOutput(outputId = "distPlot"),
        tableOutput('table')        
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot


  output$distPlot <- renderPlot({

    R = as.numeric(input$read_length)
    C = input$capture_time
    D = input$decision_time
    S = input$sequencing_speed
    X = as.numeric(input$species_abundance)/100
    
    model_eqn <- function(x) {
      return ((R + C*S) / (x*R + (1-x)*D*S + C*S))
    }

    Y = model_eqn(X)

    x <- seq(0, 1, by=0.01)
    y = lapply(x, model_eqn)
    df = data.frame(x,y)

    p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
    p + stat_function(fun = model_eqn) + xlim(0,1) +
        xlab("Abundance") +
        ylab("Enrichment") +
        annotate("Point",x=X, y=Y, colour="blue", size=3) + 
        geom_text(aes(x=X, y=Y, fontface="bold", label=paste("x=", X, "\ny=", formatC(Y, digits = 3, format = "f"))), nudge_y=0.75)
    })

  output$table <- renderTable({ 
      R = as.numeric(input$read_length)
      C = input$capture_time
      D = input$decision_time
      S = input$sequencing_speed

      model_eqn <- function(x) {
        return ((R + C*S) / (x*R + (1-x)*D*S + C*S))
      }
      abundances <- c(0.0001, 0.001, 0.01, 0.1, 0.5, as.numeric(input$species_abundance)/100)
      enrichments <- model_eqn(abundances)
      enriched_abundance <- abundances * enrichments
      df <- data.frame(abundances * 100, enrichments, enriched_abundance * 100)
      colnames(df) <- c("Starting Abundance (%)", "Predicted Enrichment", "Enriched Abundance (%)")

      df
      },  
       striped = TRUE, 
       bordered = TRUE,  
       hover = TRUE, 
       spacing = 'm',  
       align = 'l'
    )

}

shinyApp(ui = ui, server = server)