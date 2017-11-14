#' Start an Interactive castR Session
#'
#' Starts a Shiny application, in which the castR functionality can be accessed
#' interactively. The data is stored in the user's workspace and can be written
#' to a local device.
#' @param ... Parameters to be passed to \code{shiny::runApp()}.
#' @return Void.
#' @import shiny
#' @import data.table
#' @export

run_castR <- function(...) {
  
  options(shiny.maxRequestSize = 50 * 1024 ^ 2)
  
  shiny::runApp(shiny::shinyApp(
    ui = CASTLING_ui, server = CASTLING_server
  ), ...)
  
}
