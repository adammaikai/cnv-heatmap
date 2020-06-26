# CNV Analysis App

source("src/CNVPlottingFunctions.R")
load("data/CNVPlottingPrep.RData")

# Define UI for data upload app ----
ui <- fluidPage(
  title = "CNV Heatmap",
  sidebarLayout(
    sidebarPanel(
      radioButtons(inputId = 'genome', label = 'Select Genome',
                   c("GRCh37/hg19", "GRCh38/hg38", "GRCm38/mm10"), selected=character(0)),
      # Input: Select a .seg file ----
      fileInput("segFile", "Upload .seg File",
                multiple = FALSE,
                accept = c(".seg", ".txt")),
      hr(),
      # Input: Select a metadata file ----
      fileInput("metaFile", "Upload a Metadata File",
                multiple = FALSE,
                accept = c(".txt", ".csv", "tsv")),
      radioButtons(inputId = 'sep', label = 'Metadata Separator',
                   c(Tab='\t', Comma=',', Semicolon=';'),
                   '\t'),
      radioButtons("condMeta", "Annotation column", choices=character(0), selected=character(0)),
      hr(),
      actionButton(inputId = "plot_now", label = "Plot Heatmap"),
      downloadButton('downloadPlot',"Save Plot"),
      hr()
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("CNV Heatmap", id="plottab", plotOutput("Plot"))
        ),     
      HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", paste0(
          '<span data-display-if="',
          '$(&#39;html&#39;).attr(&#39;class&#39;)==&#39;shiny-busy&#39;',
          '">',
          '<i class="fa fa-spinner fa-pulse fa-fw" style="color:orange"></i>',
          '</span>'
        )),
      hr(),
      tabsetPanel(
        id = 'dataset',
        tabPanel("Segments File", id="segtab", DT::dataTableOutput("seg")),
        tabPanel("Metadata", id="metatab", DT::dataTableOutput("meta"))
      )
    )
  )
)


# Define server logic to read selected file ----
server <- function(input, output, session) {
  
  uploadseg <- eventReactive(input$segFile$datapath, {
    if(is.null(input$genome)){
      stop(safeError("Please select a reference genome!"))
    }
    # input$segFile will be NULL initially. After the user selects and uploads a file
    req(input$segFile)
    if(input$genome == "GRCh37/hg19"){
      coords <- chromCoordsHg19
    } else if(input$genome == "GRCh38/hg38"){
      coords <- chromCoordsHg38
    } else {
      coords <- chromCoordsMm10
    }
    tryCatch({
      seg <- readSeg(file.path = input$segFile$datapath, chromCoords = coords)},
      error = function(e) {stop(safeError(e))})
      list(seg=seg, coords=coords)
  })
  
  uploadmeta <- eventReactive(input$metaFile$datapath, {
    # input$segFile will be NULL initially. After the user selects and uploads a file
    req(input$metaFile)
    tryCatch({
      meta <- readMeta(input$metaFile$datapath, sep=input$sep)},
      error = function(e) {stop(safeError(e))})
    meta
  })
  
  
  output$seg <- DT::renderDataTable({
    seg <- uploadseg()$seg
    DT::datatable(seg)
  })
  
  output$meta <- DT::renderDataTable({
    meta <- uploadmeta()  
    DT::datatable(meta)
  })
  
  myPlot <- eventReactive(input$plot_now,{
    segUp <- uploadseg()
    seg <- segUp$seg
    coords <- segUp$coords
    if(!is.null(input$metaFile$datapath)){
      meta <- uploadmeta()
      column <- input$condMeta
    } else {
        meta <- NULL
        column <- NULL}
    print(coords)
    coordMat <- do.call(rbind, lapply(coords$chr,
                        function(i) {chr.i <- subset(coords, chr==i);
                        return(data.frame(row=seq(chr.i$start, chr.i$cumsum, by=100000), chr=i))}))
    cnvMat <- populateCNVMatrix(chromCoords = coords,
                                coordMat = coordMat,
                                segDf = seg)
    annos <- makeHeatmapAnnotations(cnvMat = cnvMat,
                                    chromCoords = coords,
                                    coordMat = coordMat,
                                    metadata = meta,
                                    column_anno = column)
    plotCNVHeatmap(cnvMat = cnvMat, annos = annos)
    
  })  
  
  observeEvent(input$metaFile, {
    meta <- uploadmeta()
    columns <- names(meta)[2:ncol(meta)]
    updateRadioButtons(session, "condMeta",
                             label = "Metadata Annotation",
                             choices = c("None", columns),
                             selected = columns[1])
  })

  output$Plot <-renderPlot({
    myPlot()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("CNV_heatmap_", input$condMeta, "_", gsub("-", "", Sys.Date()), sep=".pdf")
    },
    content = function (file) {
      pdf(file)
      print(myPlot())
      dev.off()
    }
  )
}

# Create Shiny app ----
shinyApp(ui, server)