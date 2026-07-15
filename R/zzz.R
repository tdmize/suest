.onLoad <- function(libname, pkgname) {
  options(marginaleffects_model_classes = unique(c(
    getOption("marginaleffects_model_classes"),
    "suest_model"
  )))

  namespace <- asNamespace("marginaleffects")
  registerS3method(
    "get_coef",
    "suest_model",
    get_coef.suest_model,
    envir = namespace
  )
  registerS3method(
    "set_coef",
    "suest_model",
    set_coef.suest_model,
    envir = namespace
  )
  registerS3method(
    "get_vcov",
    "suest_model",
    get_vcov.suest_model,
    envir = namespace
  )
  registerS3method(
    "get_predict",
    "suest_model",
    get_predict.suest_model,
    envir = namespace
  )
}
