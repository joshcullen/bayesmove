#' Dynamically explore tracks within Shiny app
#'
#' This Shiny application allows for the exploration of animal movement
#' patterns. Options are available to interactively filter the plotted tracks by
#' a selected time period of a given variable, which is then displayed on an
#' interactive map. Additionally, a data table is shown with options to filter
#' and export this table once satisfied.
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

               tabPanel("About",
                        fluidRow(
                          column(width = 8, offset = 2, h4(strong("How to use this app"))),
                          column(width = 8, offset = 2, p("This Shiny app is intended for use in exploration of animal movement patterns as well as for filtering data to meet the user's needs. At a minimum, the loaded data frame must have columns labeled ", code("id, x, y, date"), ", but can accommodate others as well. The ", code("date"), " column must be stored in ISO format (i.e., YYYY-MM-DD) of class ", code("POSIXct"), "."),
                                 br(),
                                 p("Under the 'Explore the data' header of the sidebar panel, select which track IDs you would like to visualize and the variable for which you would like to explore a time series. Click and drag on the lineplot to highlight the corresponding time range on the neighboring map. Additionally, a table of the data is printed beneath the map, which you can dynamically filter using sliders and other widgets under the 'Filter data table' header of the sidebar panel. This filtered table can then either be copied or downloaded (as a CSV file) through the buttons at the top of the table."),
                                 br(),
                                 p("Please report any issues on the", a(tags$strong("GitHub repo"), href="https://github.com/joshcullen/bayesmove/issues"), ". If interested in additional features as part of this app, please send suggestions to joshcullen10 [at] gmail [dot] com.")
                          )  #close column() for text
                        )  #close fluidRow
               ),  #close "About" tabPanel

               tabPanel("Explore data",
                        sidebarLayout(
                          sidebarPanel(
                            tags$h4(strong("Explore the data")),

                            #dropdown widget for selecting variable to viz time series
                            selectInput('var', label = 'Select a variable',
                                        choices = names(data)[!(names(data) %in%
                                                                  c("id","date"))],
                                        selected = names(data)[!(names(data) %in%
                                                                   c("id","date"))][1]),
                            #radio button to map either lines or points
                            radioButtons("radio", label = "Trajectory type",
                                         choices = c("lines", "points"),
                                         selected = "lines"),
                            br(),

                            tags$h4(strong("Filter the data")),

                            #widget/module for filtering by each column of dataframe
                            datamods::filter_data_ui("filtering", max_height = "500px")

                          ),  #close sidebar panel

                          mainPanel(
                            dygraphOutput("lineplot"),
                            leafletOutput('map'),
                            DT::dataTableOutput("tab")
                          )  #close main panel
                        )  #close sidebar layout
               )  #close "Explore data" tab panel
    )  #close navbar page
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

    # Filtering inputs for the app
    res_filter <- datamods::filter_data_server(
      id = "filtering",
      data = reactive(data),
      name = reactive("data"),
      vars = reactive(names(data)),
      widget_num = "slider",
      widget_date = "slider",
      label_na = "Missing"
    )





    ##########################################
    ### Generate output for 'Explore data' ###
    ##########################################



    # Lineplot
    output$lineplot<- renderDygraph({

      dat.ts<- res_filter$filtered() %>%
        dplyr::group_by(id) %>%
        dplyr::select(id, date, input$var) %>%
        dplyr::mutate(row = row_number()) %>%
        tidyr::pivot_wider(names_from = id, values_from = input$var) %>%
        dplyr::select(-row)

      dat.xts<- xts::xts(x = dat.ts[,-1], order.by = dat.ts$date)

      dygraph(dat.xts) %>%
        dySeries(strokeWidth = 1.5) %>%
        dyGroup(names(dat.xts), color = viridis::viridis(ncol(dat.xts))) %>%  #need to do this step for proper ordering of colors
        dyAxis("y", label = input$var, axisLabelFontSize = 16, axisLabelWidth = 75) %>%
        dyRangeSelector(dateWindow = NULL) %>%
        dyOptions(axisLineWidth = 1.5, drawGrid = FALSE) %>%
        dyHighlight(highlightCircleSize = 5,
                    highlightSeriesBackgroundAlpha = 0.2,
                    hideOnMouseOut = FALSE) %>%
        dyLegend(width = 270) %>%
        dyUnzoom() %>%
        dyCrosshair(direction = "vertical")

    })




    # Export reactive data for selected time window
    reacted.data<- eventReactive(list(res_filter$filtered(), input$lineplot_date_window), {
      req(input$lineplot_date_window)  #to prevent warning from 'if' expression below

      # define start and end times for filtering the data
      start<- strptime(input$lineplot_date_window[[1]], format = "%Y-%m-%dT%H:%M:%S",
                       tz = lubridate::tz(data$date))
      end<- strptime(input$lineplot_date_window[[2]], format = "%Y-%m-%dT%H:%M:%S",
                     tz = lubridate::tz(data$date))

      # subset dat.filt() by time window
      if (start == min(res_filter$filtered()$date) & end == max(res_filter$filtered()$date)) {
        res_filter$filtered()
      } else {
        subset = dplyr::filter(res_filter$filtered(), date >= start & date <= end)
        return(subset)
      }
    }) %>%
      debounce(millis = 500)  #add delay so map doesn't hang up


    # Map tracks for selected time window
    observeEvent(input$radio, {  #update based on whether plotting points or lines

      output$map <- renderLeaflet({

        # transform projection of coordinates to WGS84 and change from sf to data.frame
        dat.filt.sf<- sf::st_as_sf(res_filter$filtered(), coords = c("x","y"), crs = epsg) %>%
          sf::st_transform(4326) %>%
          dplyr::mutate(x = unlist(purrr::map(.data$geometry, 1)),
                        y = unlist(purrr::map(.data$geometry, 2))) %>%
          sf::st_drop_geometry()


        map1<- leaflet(data = dat.filt.sf,
                       options = leafletOptions(preferCanvas = TRUE)) %>%
          addProviderTiles(providers$Esri.OceanBasemap, group = "Ocean Basemap",
                           options = tileOptions(continuous_world = F)) %>%
          addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery",
                           options = tileOptions(continuous_world = F)) %>%
          addProviderTiles(providers$OpenStreetMap, group = "Open Street Map",
                           options = tileOptions(continuous_world = F)) %>%
          addMeasure(position = "topleft",
                     primaryLengthUnit = "kilometers",
                     primaryAreaUnit = "hectares",
                     activeColor = "#3D535D",
                     completedColor = "#7D4479") %>%
          addMiniMap(tiles = providers$Esri.OceanBasemap,
                     toggleDisplay = TRUE,
                     position = "bottomleft") %>%
          addScaleBar() %>%
          addLayersControl(baseGroups = c("World Imagery", "Ocean Basemap", "Open Street Map"),
                           overlayGroups = c("Tracks_full", "Tracks_filter"),
                           options = layersControlOptions(collapsed = TRUE))


        # add full-length tracks per ID
        if (input$radio == "lines") {  #if wanting to plot tracks as lines
          for (i in unique(dat.filt.sf$id)){
            map1 <- map1 %>%
              addPolylines(data = dat.filt.sf[dat.filt.sf$id == i,],
                           lng = ~x,
                           lat = ~y,
                           weight = 2,
                           color = "lightgrey",
                           opacity = 0.4)
          }
        } else {  #if wanting to plot tracks as points
          map1 <- map1 %>%
            addCircleMarkers(data = dat.filt.sf,
                             lng = ~x,
                             lat = ~y,
                             radius = 3,
                             fillColor = "lightgrey",
                             stroke = FALSE,
                             fillOpacity = 0.3)
        }

        map1  #print map

      })  #close renderLeaflet
    })  #close observeEvent for radio button



    # REACTIVELY UPDATE MAP
    observe({
      req(reacted.data()) # Do this if reacted.data() is not null


      # transform projection of coordinates to WGS84 and change from sf to data.frame
      dat.filt.sf<- sf::st_as_sf(res_filter$filtered(), coords = c("x","y"), crs = epsg) %>%
        sf::st_transform(4326) %>%
        dplyr::mutate(x = unlist(purrr::map(.data$geometry, 1)),
                      y = unlist(purrr::map(.data$geometry, 2))) %>%
        sf::st_drop_geometry()


      # Track w/in time window
      df<- reactive({
        df<- reacted.data()

        dat.sf<- sf::st_as_sf(df, coords = c("x","y"), crs = epsg) %>%
          sf::st_transform(4326) %>%
          dplyr::mutate(x = unlist(purrr::map(.data$geometry, 1)),
                        y = unlist(purrr::map(.data$geometry, 2))) %>%
          sf::st_drop_geometry()

        dat.sf

      })


      # First point of filtered track
      df.start.pt<- reactive({
        df<- reacted.data()

        dat.sf<- sf::st_as_sf(df, coords = c("x","y"), crs = epsg) %>%
          sf::st_transform(4326) %>%
          dplyr::mutate(x = unlist(purrr::map(.data$geometry, 1)),
                        y = unlist(purrr::map(.data$geometry, 2))) %>%
          sf::st_drop_geometry() %>%
          dplyr::group_by(id) %>%
          dplyr::slice(1)


        dat.sf

      })


      # Last point of filtered track
      df.end.pt<- reactive({
        df<- reacted.data()

        dat.sf<- sf::st_as_sf(df, coords = c("x","y"), crs = epsg) %>%
          sf::st_transform(4326) %>%
          dplyr::mutate(x = unlist(purrr::map(.data$geometry, 1)),
                        y = unlist(purrr::map(.data$geometry, 2))) %>%
          sf::st_drop_geometry() %>%
          dplyr::group_by(id) %>%
          dplyr::slice(n())

        dat.sf

      })




      # Clear old selection on map, and add new selection
      map2<- leafletProxy('map', data = df()) %>%
        clearControls() %>%
        clearShapes() %>%
        clearMarkers() %>%
        fitBounds(min(df()$x, na.rm = TRUE),
                  min(df()$y, na.rm = TRUE),
                  max(df()$x, na.rm = TRUE),
                  max(df()$y, na.rm = TRUE)) %>%
        addCircleMarkers(data = df.start.pt(),
                         lng = ~x,
                         lat = ~y,
                         fillColor = "#5EF230",
                         stroke = FALSE,
                         fillOpacity = 0.8) %>%
        addCircleMarkers(data = df.end.pt(),
                         lng = ~x,
                         lat = ~y,
                         fillColor = "red",
                         stroke = FALSE,
                         fillOpacity = 0.8) %>%
        addLegend(colors = viridis::viridis(dplyr::n_distinct(df()$id)),
                  labels = unique(df()$id),
                  opacity = 1)


      # add full-length tracks per ID
      if (input$radio == "lines") {  #tracks as lines

        for (i in 1:dplyr::n_distinct(df()$id)){
          map2 <- map2 %>%
            addPolylines(data = dat.filt.sf[dat.filt.sf$id == unique(dat.filt.sf$id)[i],],
                         lng = ~x,
                         lat = ~y,
                         weight = 2,
                         color = "lightgrey",
                         opacity = 0.4,
                         group = "Tracks_full")
        }

        # add time-filtered tracks per ID
        for (i in 1:dplyr::n_distinct(df()$id)){
          map2 <- map2 %>%
            addPolylines(data = df()[df()$id == unique(df()$id)[i],],
                         lng = ~x,
                         lat = ~y,
                         weight = 2,
                         color = viridis::viridis(dplyr::n_distinct(df()$id))[i],
                         opacity = 0.8,
                         group = "Tracks_filter")
        }

      } else {  #tracks as points

        # add full-length tracks per ID
        for (i in 1:dplyr::n_distinct(df()$id)){
          map2 <- map2 %>%
            addCircleMarkers(data = dat.filt.sf[dat.filt.sf$id == unique(dat.filt.sf$id)[i],],
                             lng = ~x,
                             lat = ~y,
                             radius = 3,
                             fillColor = "lightgrey",
                             stroke = FALSE,
                             fillOpacity = 0.3,
                             group = "Tracks_full")
        }

        # add time-filtered tracks per ID
        for (i in 1:dplyr::n_distinct(df()$id)){
          map2 <- map2 %>%
            addCircleMarkers(data = df()[df()$id == unique(df()$id)[i],],
                             lng = ~x,
                             lat = ~y,
                             radius = 3,
                             fillColor = viridis::viridis(dplyr::n_distinct(df()$id))[i],
                             stroke = FALSE,
                             fillOpacity = 0.6,
                             group = "Tracks_filter")
        }
      }

      map2  #print updated map

    })



    # Table
    output$tab<- DT::renderDataTable(
      res_filter$filtered(),   #filtered data goes here

      server = FALSE,
      extensions = "Buttons",
      options = list(paging = TRUE,
                     scrollX=TRUE,
                     searching = TRUE,
                     ordering = TRUE,
                     dom = 'Bfrtip',
                     buttons = c('copy', 'csv'),
                     pageLength=10,
                     lengthMenu=c(10,25,100) )
    )


  }
}







#################
#### Run App ####
#################

shinyApp(ui = ui(data = data, epsg = epsg),
         server = server(data = data, epsg = epsg))

}
