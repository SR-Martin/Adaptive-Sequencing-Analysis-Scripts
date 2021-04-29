library(shiny)
library(ggplot2)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Adaptive Sequencing Enrichment by Composition Model"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Slider for the number of bins ----
      textInput(inputId = "read_length", label = "Average read Length:", value = "5000"),
      sliderInput(inputId = "sequencing_speed", label = "Sequencing Speed (bps)", min = 0,  max = 1000, value = 420),
      sliderInput(inputId = "decision_time", label = "Decision Time (s)", min = 0,  max = 5, value = 1, step = 0.1),
      sliderInput(inputId = "capture_time", label = "Capture Time (s)", min = 0,  max = 5, value = 0.5, step = 0.1),


    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "distPlot")

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

    model_eqn <- function(x) {
      return ((R + C*S) / (x*R + (1-x)*D*S + C*S))
    }

    x <- seq(0, 1, by=0.01)
    y = lapply(x, model_eqn)
    df = data.frame(x,y)

    p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
    p + stat_function(fun = model_eqn) + xlim(0,1) +
        labs(title = "Predicted Enrichment") +
        xlab("Abundance") +
        ylab("Enrichment")
    })

}

shinyApp(ui = ui, server = server)