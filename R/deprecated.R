## Deprecated names ----------------------------------------------------------
##
## Renamed to follow the package-wide camelCase convention. The old names keep
## working but signal via .Deprecated() and forward to the replacement. Each
## alias reproduces the *old* signature and forwards by name, so that neither
## positional nor named calls can be silently remapped onto a different formal.

#' Deprecated functions in dMod
#'
#' @description
#' These functions have been renamed to follow the package-wide camelCase
#' convention. The old names still work but signal a deprecation warning and
#' will be removed in a future release.
#'
#' \itemize{
#'   \item `distributed_computing()` -> [distributedComputing()]
#'   \item `profile_pars_per_node()` -> [profileParsPerNode()]
#'   \item `mname()` -> [modelname()]
#'   \item `import_sbml()` -> [importSbml()]
#'   \item `export_sbml()` -> [exportSbml()]
#'   \item `read_petab_yaml()` -> [readPetabYaml()]
#'   \item `read_petab_tables()` -> [readPetabTables()]
#' }
#'
#' @name dMod-deprecated
#' @keywords internal
NULL


#' @rdname dMod-deprecated
#' @inheritParams distributedComputing
#' @export
distributed_computing <- function(...,
                                  jobname,
                                  partition = "single",
                                  cores = 16,
                                  nodes = 1,
                                  mem_per_core = 2,
                                  walltime = "01:00:00",
                                  ssh_passwd = NULL,
                                  machine = "cluster",
                                  var_values = NULL,
                                  no_rep = NULL,
                                  recover = TRUE,
                                  purge_local = FALSE,
                                  compile = FALSE,
                                  link = FALSE,
                                  custom_folders = NULL,
                                  resetSeeds = TRUE,
                                  returnAll = TRUE) {
  .Deprecated("distributedComputing", package = "dMod")
  ## The `...` of distributedComputing() is captured *unevaluated* via
  ## substitute(); forwarding dots preserves the original promise expressions,
  ## so the remote code arrives intact.
  distributedComputing(...,
                       jobname = jobname,
                       partition = partition,
                       cores = cores,
                       nodes = nodes,
                       mem_per_core = mem_per_core,
                       walltime = walltime,
                       ssh_passwd = ssh_passwd,
                       machine = machine,
                       var_values = var_values,
                       no_rep = no_rep,
                       recover = recover,
                       purge_local = purge_local,
                       compile = compile,
                       link = link,
                       custom_folders = custom_folders,
                       resetSeeds = resetSeeds,
                       returnAll = returnAll)
}


#' @rdname dMod-deprecated
#' @inheritParams profileParsPerNode
#' @export
profile_pars_per_node <- function(parameters, fits_per_node,
                                  side = c("both", "split")[1]) {
  .Deprecated("profileParsPerNode", package = "dMod")
  profileParsPerNode(parameters = parameters,
                     fits_per_node = fits_per_node,
                     side = side)
}


#' @rdname dMod-deprecated
#' @param x dMod object
#' @param conditions character vector of conditions
#' @export
mname <- function(x, conditions = NULL) {
  .Deprecated("modelname", package = "dMod")
  modelname(x, conditions = conditions)
}


#' @rdname dMod-deprecated
#' @inheritParams importSbml
#' @export
import_sbml <- function(modelpath) {
  .Deprecated("importSbml", package = "dMod")
  importSbml(modelpath = modelpath)
}


#' @rdname dMod-deprecated
#' @inheritParams exportSbml
#' @export
export_sbml <- function(eqnlist, parameters = NULL, inits = NULL, filepath,
                        modelID = "dMod_export") {
  .Deprecated("exportSbml", package = "dMod")
  exportSbml(eqnlist = eqnlist, parameters = parameters, inits = inits,
             filepath = filepath, modelID = modelID)
}


#' @rdname dMod-deprecated
#' @inheritParams readPetabYaml
#' @export
read_petab_yaml <- function(yamlPath) {
  .Deprecated("readPetabYaml", package = "dMod")
  readPetabYaml(yamlPath = yamlPath)
}


#' @rdname dMod-deprecated
#' @inheritParams readPetabTables
#' @export
read_petab_tables <- function(yamlPath) {
  .Deprecated("readPetabTables", package = "dMod")
  readPetabTables(yamlPath = yamlPath)
}
