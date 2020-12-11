#' Dynamically explore tracks within Shiny app
#'
#' This Shiny application allows for the exploration of individual movement
#' patterns, as well as those of all tracked individuals. Options are available
#' to interactively filter the plotted tracks by a selected time period of a
#' given variable. Additionally, collections of individuals can be plotted
#' simultaneously and these tracks can be filtered to only display movements
#' within a given period of time.
#'
#' Currently, the time series plot shown for the exploration of individual
#' tracks cannot display variables of class \code{character} or \code{factor}.
#' Therefore, these should be changed to numeric values if they are to be
#' plotted.
#'
#' If the data are stored as longitude and latitude (i.e., WGS84), the EPSG code
#' is 4326. All other codes will need to be looked up if they are not already
#' known.
#'
#' @param data A data frame that must contain columns labeled \code{id, x, y,
#'   date}, but can include any other variables of interest.
#' @param epsg numeric. The coordinate reference system (CRS) as an EPSG code.
#'
#' @examples
#' \dontrun{
#' #load data
#' data(tracks)
#'
#' #run Shiny app
#' shiny_tracks(data = tracks, epsg = 32617)
#'
#' }
#'
#' @import shiny
#' @import dygraphs
#' @import leaflet
#'
#'
#' @export
shiny_tracks = function(data, epsg){



###################
#### Define UI ####
###################

ui <- function(data, epsg) {
  navbarPage("Exploration of Animal Telemetry Data",
             theme = shinythemes::shinytheme("flatly"),
             tabPanel("Explore data",
                      sidebarLayout(
                        sidebarPanel(selectInput('animal_id', label = 'Select an ID',
                                                 choices = unique(data$id),
                                                 selected = unique(data$id)[1]),
                                     selectInput('var', label = 'Select a Variable',
                                                 choices = names(data)[names(data) != "id"],
                                                 selected = names(data)[names(data) != "id"][1]),
                                     p("Use the ", strong("Explore data"), " tab to explore temporal patterns of user selected variables from the loaded data frame for each ID in the dataset. At minimum, the data frame must have columns labeled ", code("id, x, y, date"), ", but can accommodate others as well. The ", code("date"), " column must be stored in ISO format (i.e., YYYY-MM-DD) of class ", code("POSIXct"), ". Click and drag on the lineplot to highlight the corresponding time range on the neighboring map.", style = "font-family: 'times'; font-si16pt"),
                                     p("Click the ", strong("View all tracks"), " tab to view tracks for all individuals. This provides a quick overview of all track segments separately and allows for filtering tracks by time.", style = "font-family: 'times'; font-si16pt"),
                                     br(),
                                     br(),
                                     p("Application author: ",
                                       a("Josh Cullen", href = "https://joshcullen.github.io/"),

                                       tags$br(),

                                       a("University of Florida", href = "http://www.ufl.edu"),
                                       style = "font-family: 'times'; font-si16pt")
                        ),

                        mainPanel(dygraphOutput("lineplot"),
                                  leafletOutput('map'))
                      )
             ),
             tabPanel("View all tracks",
                      sidebarLayout(
                        sidebarPanel(
                          shinyWidgets::pickerInput(
                            inputId = "select_ids",
                            label = "Select/deselect by ID",
                            choices = unique(data$id),
                            selected = unique(data$id),
                            options = list(
                              `actions-box` = TRUE),
                            multiple = TRUE
                          ),
                          p("Use the ", strong("Explore data"), " tab to explore temporal patterns of user selected variables from the loaded data frame for each ID in the dataset. At minimum, the data frame must have columns labeled ", code("id, x, y, date"), ", but can accommodate others as well. The ", code("date"), " column must be stored in ISO format (i.e., YYYY-MM-DD) of class ", code("POSIXct"), ". Click and drag on the lineplot to highlight the corresponding time range on the neighboring map.", style = "font-family: 'times'; font-si16pt"),
                          p("Click the ", strong("View all tracks"), " tab to view tracks for all individuals. This provides a quick overview of all track segments separately and allows for filtering tracks by time.", style = "font-family: 'times'; font-si16pt"),
                          br(),
                          br(),
                          p("Application author: ",
                            a("Josh Cullen", href = "https://joshcullen.github.io/"),

                            tags$br(),

                            a("University of Florida", href = "http://www.ufl.edu"),
                            style = "font-family: 'times'; font-si16pt")),

                        mainPanel(sliderInput(inputId = "range", label = NULL,
                                                  min = min(data$date),
                                                  max = max(data$date),
                                                  value = range(data$date), width = '100%',
                                                  step = 90),
                                  leafletOutput('map_all'))
                      )
             )
  )
}





#######################
#### Define Server ####
#######################

server <- function(data, epsg) {
  function(input, output, session) {

    #####################
    #### Modify data ####
    #####################

    data$id<- as.character(data$id)
    data$date<- lubridate::as_datetime(data$date)


    #############################################
    ### Update sidebar and make reactive data ###
    #############################################

    ### 'Explore data' tab

    # Filter by ID
    dat.filt <- reactive({
      d<- data[data$id == input$animal_id, ]

      return(d)
    })




    ### 'View all tracks' tab

    # Make reactive data for selected IDs
    dat.list<- reactive({
      d.list<- split(data, data$id)  #split data into list by ID
      d2<- dplyr::bind_rows(d.list[input$select_ids])  #only keep selected IDs

      return(d2)
    })

    #make reactive data based on slider input
    dat.filt2<- reactive({
      dat<- dat.list()
      dat<- dat[dat$date >= input$range[1] & dat$date <= input$range[2], ]  #only keep w/in range

      #format data where each ID is a linestring
      dat.sf<- sf::st_as_sf(dat, coords = c("x","y"), crs = epsg) %>%
        sf::st_transform(4326) %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(do_union=F) %>%
        sf::st_cast("LINESTRING")

      return(dat.sf)
    })






    ##########################################
    ### Generate output for 'Explore data' ###
    ##########################################

    # Lineplot
    output$lineplot<- renderDygraph({

      dygraph(xts::xts(x = dat.filt()[,input$var], order.by = dat.filt()$date)) %>%
        dySeries(label = input$var, strokeWidth = 1.5) %>%
        dyAxis("y", label = input$var, axisLabelFontSize = 16, axisLabelWidth = 75) %>%
        dyRangeSelector(dateWindow = NULL) %>%
        dyOptions(axisLineWidth = 1.5, drawGrid = FALSE, colors = "black") %>%
        dyLegend(width = 270) %>%
        dyUnzoom() %>%
        dyCrosshair(direction = "vertical")

    })


    # Export reactive data for selected time window
    reacted.data<- reactive({
      req(input$lineplot_date_window)  #to prevent warning from 'if' expression

      start=strptime(input$lineplot_date_window[[1]], format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
      end=strptime(input$lineplot_date_window[[2]], format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")

      if (start == min(dat.filt()$date) & end == max(dat.filt()$date)){
        dat.filt()
      } else {
        subset = dplyr::filter(dat.filt(), date >= start & date <= end)
        return(subset)
      }
    })

    # Map tracks for selected time window
    output$map <- renderLeaflet({

      dat.filt.sf<- sf::st_as_sf(dat.filt(), coords = c("x","y"), crs = epsg) %>%
        sf::st_transform(4326) %>%
        sf::st_cast("LINESTRING")


      leaflet(data = dat.filt.sf, options = leafletOptions(preferCanvas = TRUE)) %>%
        addProviderTiles(providers$Esri.OceanBasemap, group = "Ocean Basemap",
                         options = tileOptions(continuous_world = F)) %>%
        addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery",
                         options = tileOptions(continuous_world = F)) %>%
        addProviderTiles(providers$CartoDB.DarkMatterNoLabels, group = "Dark Map",
                         options = tileOptions(continuous_world = F)) %>%
        addPolylines(lng = as.numeric(sf::st_coordinates(dat.filt.sf)[,1]),
                     lat = as.numeric(sf::st_coordinates(dat.filt.sf)[,2]),
                     weight = 2,
                     color = "lightgrey",
                     opacity = 0.4) %>%
        addMeasure(position = "topleft",
                   primaryLengthUnit = "kilometers",
                   primaryAreaUnit = "hectares",
                   activeColor = "#3D535D",
                   completedColor = "#7D4479") %>%
        addMiniMap(tiles = providers$Esri.OceanBasemap,
                   toggleDisplay = TRUE,
                   position = "bottomleft") %>%
        addScaleBar() %>%
        addLayersControl(baseGroups = c("World Imagery", "Ocean Basemap", "Dark Map"),
                         options = layersControlOptions(collapsed = TRUE))
    })




    observe({

      # UPDATED MAP
      req(reacted.data()) # Do this if reacted.data() is not null


      dat.filt.sf<- sf::st_as_sf(dat.filt(), coords = c("x","y"), crs = epsg) %>%
        sf::st_transform(4326) %>%
        sf::st_cast("LINESTRING")

      # Track w/in time window
      df<- reactive({
        df<- reacted.data()

        dat.sf<- sf::st_as_sf(df, coords = c("x","y"), crs = epsg) %>%
          sf::st_transform(4326) %>%
          sf::st_cast("LINESTRING")

        dat.sf

      })

      # First point of filtered track
      df.start.pt<- reactive({
        df<- reacted.data()

        dat.sf<- sf::st_as_sf(df, coords = c("x","y"), crs = epsg) %>%
          sf::st_transform(4326) %>%
          dplyr::slice(1)

        dat.sf

      })

      # Last point of filtered track
      df.end.pt<- reactive({
        df<- reacted.data()

        dat.sf<- sf::st_as_sf(df, coords = c("x","y"), crs = epsg) %>%
          sf::st_transform(4326) %>%
          dplyr::slice(nrow(df))

        dat.sf

      })




      # Clear old selection on map, and add new selection
      leafletProxy('map', data = df()) %>%
        clearControls() %>%
        clearShapes() %>%
        clearMarkers() %>%
        fitBounds(as.numeric(sf::st_bbox(df())[1]),
                  as.numeric(sf::st_bbox(df())[2]),
                  as.numeric(sf::st_bbox(df())[3]),
                  as.numeric(sf::st_bbox(df())[4])) %>%
        addPolylines(lng = as.numeric(sf::st_coordinates(dat.filt.sf)[,1]),
                     lat = as.numeric(sf::st_coordinates(dat.filt.sf)[,2]),
                     weight = 2,
                     color = "lightgrey",
                     opacity = 0.4) %>%
        addPolylines(lng = as.numeric(sf::st_coordinates(df())[,1]),
                     lat = as.numeric(sf::st_coordinates(df())[,2]),
                     weight = 2,
                     color = "darkturquoise",
                     opacity = 0.8) %>%
        addCircleMarkers(data = df.start.pt(),
                         fillColor = "#5EF230",
                         stroke = FALSE,
                         fillOpacity = 0.8) %>%
        addCircleMarkers(data = df.end.pt(),
                         fillColor = "red",
                         stroke = FALSE,
                         fillOpacity = 0.8)

    })








    #############################################
    ### Generate output for 'View all tracks' ###
    #############################################

    # Create map
    output$map_all <- renderLeaflet({

      #set palette
      pal1<- colorFactor(palette = "viridis", domain = dat.filt2()$id)

      labs<- sprintf("<strong>ID %s</strong><br/>", dat.filt2()$id) %>%
        lapply(htmltools::HTML)


      leaflet(data = dat.filt2()) %>%
        addProviderTiles(providers$Esri.OceanBasemap, group = "Ocean Basemap",
                         options = tileOptions(continuous_world = F)) %>%
        addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery",
                         options = tileOptions(continuous_world = F)) %>%
        addProviderTiles(providers$CartoDB.DarkMatterNoLabels, group = "Dark Map",
                         options = tileOptions(continuous_world = F)) %>%
        clearControls() %>%
        clearShapes() %>%
        fitBounds(as.numeric(sf::st_bbox(dat.filt2())[1]),
                  as.numeric(sf::st_bbox(dat.filt2())[2]),
                  as.numeric(sf::st_bbox(dat.filt2())[3]),
                  as.numeric(sf::st_bbox(dat.filt2())[4])) %>%
        addPolylines(data = dat.filt2(),
                     weight = 2,
                     color = ~pal1(id),
                     label = labs) %>%
        leaflet::addLegend("bottomright",
                           pal = pal1,
                           values = dat.filt2()$id,
                           title = "ID",
                           opacity = 1) %>%
        addMiniMap(tiles = providers$Esri.OceanBasemap,
                   toggleDisplay = TRUE,
                   position = "bottomleft") %>%
        addScaleBar() %>%
        addLayersControl(baseGroups = c("Ocean Basemap", "World Imagery", "Dark Map"),
                         options = layersControlOptions(collapsed = TRUE), position = "topleft")
    })

  }
}







#################
#### Run App ####
#################

shinyApp(ui = ui(data = data, epsg = epsg),
         server = server(data = data, epsg = epsg))

}
