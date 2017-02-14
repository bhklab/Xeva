##------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------
#' Fit a GP using Python's GPy package.
#'
#' \code{fitGP} fit a GP using Python's GPy package
#'
#' @param replicates A vector of multiple replicates that user would like to fit GP on
#' @param source Data source
#'
#' @return results A list of time, predicted mean, and variance of fit GP
#'
#' @examples
#' fitGP(c("PHLC111_P7.701.A1", "PHLC111_P7.703.A3", "PHLC111_P7.706.B1", "PHLC111_P7.708.B3"), lpdx)
#' fitGP(c("PHLC111_P7.702.A2", "PHLC111_P7.704.A4", "PHLC111_P7.705.A5", "PHLC111_P7.707.B2"), lpdx)
#' @export
fitGP <- function(replicates, source) {
  times <- c()
  volumes <- c()

  max_time_index <- 0
  max_volume_index <- 0

  for (i in 1:length(replicates)) {
    df <- getExperiment(source, replicates[i])
    times[[i]] <- df$time
    volumes[[i]] <- df$volume

    if (length(df$time) > max_time_index) {
      max_time_index <- length(df$time)
    }

    if (length(df$volume) > max_volume_index) {
      max_volume_index <- length(df$volume)
    }
  }

  # reformat times and volumes
  for (i in 1:length(times)) {
    times[[i]] <- c(times[[i]], rep(NA, max_time_index - length(times[[i]])))
    volumes[[i]] <- c(volumes[[i]], rep(NA, max_volume_index - length(volumes[[i]])))
  }

  # write all to tempfiles
  write.table(data.frame(lapply(t(times), as.numeric)),
              file="python/temptime.csv", row.names = FALSE, col.names = FALSE)
  write.table(data.frame(lapply(volumes, as.numeric), stringsAsFactors = FALSE),
              file="python/tempvolume.csv", row.names = FALSE, col.names = FALSE)

  # run python script
  system("python3 python/fit_single_gp.py")

  # get results
  results = read.csv("python/tempresults.csv", sep = " ")

  return(results)
}

##------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------
#' Get KL divergence of two GPs
#'
#' \code{getGPStatistics} Get the KL divergence of two GPs.
#'
#' @param controlReplicates Vector of control replicates
#' @param caseReplicates Vector of case replicates
#' @param source Data source of replicates
#'
#' @return  Returns the KL divergence value for two GPs.
#'
#' @examples
#' getGPStatistics(c("PHLC111_P7.701.A1", "PHLC111_P7.703.A3"), c("PHLC111_P7.705.A5", "PHLC111_P7.707.B2"), lpdx)
#'
#' @export
getGPStatistics <- function(controlReplicates, caseReplicates, source) {
  control_times <- c()
  case_times <- c()

  control_volumes <- c()
  case_volumes <- c()

  max_time_index <- 0
  max_volume_index <- 0

  for (i in 1:length(controlReplicates)) {
    df <- getExperiment(source, controlReplicates[i])
    control_times[[i]] <- df$time
    control_volumes[[i]] <- df$volume

    if (length(df$time) > max_time_index) {
      max_time_index <- length(df$time)
    }

    if (length(df$volume) > max_volume_index) {
      max_volume_index <- length(df$volume)
    }
  }

  for (i in 1:length(caseReplicates)) {
    df <- getExperiment(source, caseReplicates[i])
    case_times[[i]] <- df$time
    case_volumes[[i]] <- df$volume

    if (length(df$time) > max_time_index) {
      max_time_index <- length(df$time)
    }

    if (length(df$volume) > max_volume_index) {
      max_volume_index <- length(df$volume)
    }
  }

  # reformat times and volumes
  for (i in 1:length(control_times)) {
    control_times[[i]] <- c(control_times[[i]], rep(NA, max_time_index - length(control_times[[i]])))
    control_volumes[[i]] <- c(control_volumes[[i]], rep(NA, max_volume_index - length(control_volumes[[i]])))
  }

  # reformat times and volumes
  for (i in 1:length(case_times)) {
    case_times[[i]] <- c(case_times[[i]], rep(NA, max_time_index - length(case_times[[i]])))
    case_volumes[[i]] <- c(case_volumes[[i]], rep(NA, max_volume_index - length(case_volumes[[i]])))
  }

  # write all to tempfiles
  write.table(data.frame(lapply(t(control_times), as.numeric)),
              file="python/temp_control_time.csv", row.names = FALSE, col.names = FALSE)
  write.table(data.frame(lapply(control_volumes, as.numeric), stringsAsFactors = FALSE),
              file="python/temp_control_volume.csv", row.names = FALSE, col.names = FALSE)
  write.table(data.frame(lapply(t(case_times), as.numeric)),
              file="python/temp_case_time.csv", row.names = FALSE, col.names = FALSE)
  write.table(data.frame(lapply(case_volumes, as.numeric), stringsAsFactors = FALSE),
              file="python/temp_case_volume.csv", row.names = FALSE, col.names = FALSE)

  # run python script
  system("python3 python/two_gp_statistics.py")

  # get results
  results = read.csv("python/temp_statistics_results.csv", sep = " ")

  return(results)
}
