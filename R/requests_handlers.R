#' Repeat failing functions and manage errors.
#'
#' @param fun The same function supplied to ncbi_limit_handler.
#' @param attempts The number of times to repeat the function before stopping the while loop.
#'
#' @return The result of the function provided.
#' @export
#'
#' @keywords internal
#'
#' @description
#' This function, although used internally, is exported in order to be able to
#' be called by `ncbi_limit_handler` inside its future::future function.
#'
connection_handler <- function(fun, attempts = 5) {

  lock <- 0

  while (attempts > 0) {

    result <- tryCatch({

      withCallingHandlers({

        fun()

      }, warning = function(w) {

        lock <<- 1

        warn <- paste("WARNING:", conditionMessage(w))

        handled <- attempts_handler(warn, attempts)

        attempts <<- handled[[1]]

        if (attempts == 0) {
          lock <<- 0
        }

        if (!is.null(handled[[2]])) {
          stop(handled[[2]])
        }

      })

    }, error = function(e){

      paste("ERROR:", conditionMessage(e))

    })

    if (lock == 0) {

      handled <- attempts_handler(result, attempts)

      attempts <- handled[[1]]

    }

  }

  result

}

#' Run function on a list/vector of elements following a specific timing in a future framework.
#'
#' @param data The list/vector whose elements will be provided to the function in the 'fun' parameter.
#' @param message the message to print with the progress bar.
#' @param api_rate The rate with which to perform the function on each element. Must be a number between 3 and 10 which will translate in a rate of '1 / api_rate' seconds.
#' @param fun The function to run on each element of the list/vector at the 'elements' parameter.
#' @param seed due to the functioning of the bold package, this allows to ignore random number generation warnings by the future package. Keep to FALSE (default of future::future function) for all other cases
#'
#' @return The result of the function provided with the fun parameter.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' This function allows to run another function as a future (future::future()) object. If the current future plan is "sequential" (as in every 'normal' condition of NON-parallelization), each element of the list/vector at the elements parameter.
#' If, on the other hand, the future plan has been previously set to "multisession", then the function applied to each element will be sent in background as a future object stored in a list, allowing to apply the function to multiple elements each 1 / api_rate seconds. When all elements have been processed with the fun, then future::value will extract the results and return them as a list.
#' This allows each user to adopt the wanted future plan, eventually allowing to speed up the analyses with a dynamic parallelization that respects the http requests limit imposed by the NCBI.
#' This function also implements the other function connection_handler which allows to repeat the function if there is a connection error.
#'
#' @examples \dontrun{
#' ids <- c("Polymastia invaginata", "Aglaophamus trissophyllus", "Hastingsia gracilis")
#'
#' ncbi_limit_handler(ids, api_rate = 3, function(id){
#' ncbi_searcher(id, "taxonomy")
#' })
#' }
ncbi_limit_handler <- function(data, api_rate = 3, fun, message = "LoremIpsum", seed = FALSE) {

  p <- progressr::progressor(steps = length(data))

  elements <- 1:length(data)

  results_list <- list()
  pos <- 0

  for (el in elements) {

    connection_handler
    attempts_handler

    Sys.sleep(1 / api_rate)

    pos <- pos + 1
    args <- as.list(stats::setNames(el, nm = methods::formalArgs(fun)))

    if (api_rate > 3) {

      key <- Sys.getenv("ENTREZ_KEY")

    }

    results_list[[pos]] <- future::future({

      if (api_rate > 3) {

        Sys.setenv(ENTREZ_KEY = key)

      }

      result <- connection_handler(function() {
        return(do.call(fun, args))
      })

      if (length(grep("ERROR:", result)) > 0) {

        values <- data[el]

        stop(result, "\nOccurred at position '", el, "'\nWith values:\n", values)

      } else {

        result

      }

    }, earlySignal = TRUE, seed = seed)

    if (!is.null(p)) {

      p(message = sprintf(message))

    }

  }

  future::resolve(results_list)

  future::value(results_list)

}

#' Define errors and warnings and choose proper repeating scheme accordingly.
#'
#' @param input the warning/error as captured by tryCatch and withCallingHandlers
#' @param attempts the current number of attempts as computed by connection_handler
#'
#' @return a number defining the correct attempts value
#'
#' @keywords internal
#' @noRd
#'
attempts_handler <- function(input, attempts) {

  # if no errors or warnings have been produced, change attempts to 0 in order
  # to exit the while loop
  if (length(grep("ERROR:", input)) == 0 &&
      length(grep("WARNING:", input)) == 0) {

    return(list(0, NULL))

  }

  # if the error/warning contains this specific message, wait a bit more and
  # repeat 10 times, then exit and print error/warning
  if (length(grep("Can't fetch uids from history because of: Result not ready",
                  input)) > 0) {

    Sys.sleep(1)

    if (attempts < 10) {

      attempts <- attempts + 1

      return(list(attempts, NULL))

    } else {

      return(list(0, NULL))

    }

  }

  # if the error/warning contains this specific message, do not repeat, as it
  # won't affect the outcome
  if (length(grep("In file.exists(file) : unable to translate",
                  input)) > 0) {

    return(list(0, NULL))

  }

  # this warning comes from bold_tax_name searches and indicates the bold servers
  # have blocked the user. In this case, turn the warning into an error and stop
  if (length(grep("Content was type '' when it should've been type 'text/html; charset=utf-8'",
                  input)) > 0) {

    return(list(0, "Blocked by BOLD servers"))

  }

  # similarly to the above error, this indicates blocking by BOLD servers
  if (length(grep("You have exceeded your allowed",
                  input)) > 0) {

    return(list(0, "Blocked by BOLD servers"))

  }

  # similarly to the above error, this indicates blocking by BOLD servers
  if (length(grep("DOCTYPE html",
                  input)) > 0) {

    return(list(0, "Blocked by BOLD servers"))

  }

  # sometimes the NCBI blocks temporarily
  if (length(grep("Could not resolve host",
                  input)) > 0) {

    return(list(0, "NCBI servers unavailable"))

  }

  # for other generic errors/warnings, wait one second and remove one attempt
  if ((length(grep("ERROR:", input)) > 0 ||
       length(grep("WARNING:", input)) > 0)) {

    Sys.sleep(1)

    attempts <- attempts - 1

    return(list(attempts, NULL))

  }

}
