#' Customized Factory for Learner Wrappers
#'
#' This function is a modified version of `create.Learner` from the
#' SuperLearner package (Polley et al., 2011).
#' @references
#' Polley et al. (2011).
#' SuperLearner: Super Learner Prediction. URL: \url{https://CRAN.R-project.org/package=SuperLearner}
#'
#' @export
create.Learner = function(base_learner, params = list(), tune = list(),
                                   env = parent.frame(), name_prefix = base_learner,
                                   detailed_names = F, verbose = F) {
  if (length(tune) > 0) {
    tuneGrid = expand.grid(tune, stringsAsFactors = FALSE)
    nameGrid = tuneGrid
    if ("num_intervals" %in% colnames(nameGrid)) {
      nameGrid$num_intervals <- paste0("B", nameGrid$num_intervals)
    }
    names = rep("", nrow(nameGrid))
    max_runs = nrow(tuneGrid)
  } else {
    # Run once if no tuneGrid is defined, otherwise run once per grid row.
    max_runs = 1
    tuneGrid = NULL
    names = c()
  }

  for (i in seq(max_runs)) {
    name = paste(name_prefix, i, sep="_")

    if (length(tune) > 0) {
      # Specify drop=F in case tuneGrid is a single-column dataframe.
      g = tuneGrid[i, , drop=F]
      h = nameGrid[i, , drop=F]
      g = format(g, scientific = FALSE)
      h = format(h, scientific = FALSE)
      # Separate with "_" because some hyperparameters could be floats with a period.
      if (detailed_names) {
        name = do.call(paste, c(list(name_prefix), h, list(sep="_")))
      }
    } else {
      g = c()
    }

    names[i] = name

    # Create the custom learner function. This approach allows us to not specify
    # some of the learner arguments. and have the function use its own defaults.
    # Or we can set those arguments to "NULL".
    fn_params = ""
    all_params = c(as.list(g), params)
    for (name_i in names(all_params)) {
      val = all_params[[name_i]]
      if (!is.null(val) && all(val != "NULL")) {
        if (is.vector(val) && !is.character(val)) {
          val <- paste0("c(", paste(val, collapse = ", "), ")")
        } else if (is(val, "character")) {
          val <- paste0('"', val, '"')
        }
        fn_params <- paste0(fn_params, ", ", name_i, "=", val)
      }
    }
    fn = paste0(name, " <- function(...) ", base_learner, "(...", fn_params, ")")
    if (verbose) {
      cat(fn, "\n")
    }
    eval(parse(text = fn), envir = env)
  }
  results = list(grid = tuneGrid, names = names, base_learner = base_learner,
                 params = params)
  invisible(results)
}
