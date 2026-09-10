dat <- read.csv("suest_r_panel_poisson_re_crosslang_benchmark.csv")

lower_matrix <- function(values, size) {
  out <- matrix(0, size, size)
  position <- 1L
  for (row in seq_len(size)) {
    for (column in seq_len(row)) {
      out[row, column] <- values[position]
      position <- position + 1L
    }
  }
  out + t(out) - diag(diag(out))
}

model_reference <- function(coefficients, covariance, logLik) {
  order <- c(3L, 1L, 2L, 4L)
  parameter_names <- c("(Intercept)", "x", "z", "lnalpha")
  coefficients <- coefficients[order]
  covariance <- covariance[order, order]
  names(coefficients) <- parameter_names
  dimnames(covariance) <- list(parameter_names, parameter_names)
  list(coefficients = coefficients, vcov = covariance, logLik = logLik)
}

system_reference <- function(coefficients, covariance, nobs_union,
                             nobs_overlap, n_clusters) {
  order <- c(3L, 1L, 2L, 4L, 7L, 5L, 6L, 8L)
  parameter_names <- unlist(lapply(c("Y1", "Y2"), function(model)
    paste0(model, "::", c("(Intercept)", "x", "z", "lnalpha"))))
  coefficients <- coefficients[order]
  covariance <- covariance[order, order]
  names(coefficients) <- parameter_names
  dimnames(covariance) <- list(parameter_names, parameter_names)
  list(
    coefficients = coefficients,
    vcov = covariance,
    nobs_union = nobs_union,
    nobs_overlap = nobs_overlap,
    n_clusters = n_clusters
  )
}

balanced_model1 <- model_reference(
  c(0.364062082090571, -0.201215659844019, 0.231569576644006,
    -0.820824998084864),
  lower_matrix(c(
    0.00184357772607491,
    -0.000126193370586697, 0.00186301875514393,
    -0.000595665049902791, 0.000266607573038125, 0.00728158339947888,
    0.00028223538734099, -0.000302724614993045, -0.0000129661039540519,
      0.0439275450180362
  ), 4L),
  -705.545047704938
)

balanced_model2 <- model_reference(
  c(0.244058623712813, 0.116956410427222, -0.142242745897189,
    -0.862056924065249),
  lower_matrix(c(
    0.00248636253281761,
    -0.0000864782997740059, 0.00261429365342762,
    -0.000729386754403628, -0.000372843521719395, 0.0078339532043387,
    -0.00074793177798095, -0.00072090188563227, 0.000115211581476923,
      0.0530519476974698
  ), 4L),
  -583.866264482637
)

balanced_system <- system_reference(
  c(0.364062082090571, -0.201215659844019, 0.231569576644006,
    -0.820824998084864, 0.244058623712813, 0.116956410427222,
    -0.142242745897189, -0.862056924065249),
  lower_matrix(c(
    0.00202692088466012,
    -0.000490136695006383, 0.00267709070873796,
    -0.00136429949007211, 0.000886645348499617, 0.00726673095225591,
    0.0124216827778051, -0.0103246501764487, -0.0179514141410127,
      1.41396833812739,
    -0.000263697932367644, -0.0000766446824354066,
      -0.0000127813035694047, 0.00854259219487706, 0.0025061264266764,
    0.000150785200976754, -0.0000936403116653143,
      -0.000584836045282156, 0.00280406803186405,
      0.000298059378431292, 0.00273170642914627,
    0.0000682382543519136, 0.000296166014508419,
      0.00242893797166897, 0.00327621176609108, -0.00203110251804999,
      -0.000996871848516599, 0.00940192750106368,
    0.0100961563477326, 0.00581298803228583, 0.00210764974142037,
      0.0569513361413991, -0.0235132156929146, -0.0258877646236273,
      0.031993569869325, 1.87314517015844
  ), 8L),
  480L, 480L, 80L
)

partial_system <- system_reference(
  c(0.370590300877418, -0.199442278067673, 0.229122855932255,
    -0.842609683645375, 0.227548452125264, 0.154346557657659,
    -0.163450275493318, -0.921293783290415),
  lower_matrix(c(
    0.00696384877697222,
    -0.00118609271959331, 0.00422405705817986,
    -0.00219780412611074, 0.00100948688408879, 0.00107514187399745,
    0.0537527477294236, -0.0459106386205527, 0.00392086236044623,
      5.99983144677677,
    0.000906078973152929, 0.000537488518539954,
      -0.000410870403820952, -0.0431782221336795, 0.00361515330146956,
    -0.0000836556510398445, 0.00139042817282939,
      0.00010708243383759, -0.0832361823840758,
      -0.000366583207490882, 0.00332950520726374,
    -0.0000830576618205684, -0.000406661376742041,
      0.000150806674424317, 0.0123937046334123,
      -0.000757159594579274, 0.000379587265838371,
      0.000700060842622785,
    0.0511268308919713, -0.0461299212615735, -0.00164077077945218,
      5.32402334640773, -0.0315733921634079, -0.103150340044152,
      -0.0102581053164703, 6.78514875843941
  ), 8L),
  480L, 444L, 16L
)

disjoint_system <- system_reference(
  c(0.431987390204632, -0.262263382330325, 0.0674156067318728,
    -1.08920269039853, 0.311008924203921, 0.150568737069708,
    -0.137799055941583, -0.987172069117551),
  lower_matrix(c(
    0.00417460914710092,
    -0.000567315772557421, 0.00618653752776584,
    -0.00428508614041567, 0.00163008458828885, 0.0140533025584256,
    0.0349931124024754, -0.0526358110805416, -0.0494750726558251,
      2.86007936416092,
    0, 0, 0, 0, 0.00537289848197433,
    0, 0, 0, 0, -0.000717685953679169, 0.00441217467588528,
    0, 0, 0, 0, -0.00376696028325011, -0.0015161480816032,
      0.0172384674911182,
    0, 0, 0, 0, -0.0254964275741972, -0.0128615130535468,
      0.0626078536573528, 4.13795104897097
  ), 8L),
  480L, 0L, 80L
)

reference <- list(
  metadata = list(
    Stata = "19.5",
    suest2 = "1.0.0",
    benchmark_revision = 1L,
    ancillary_score_convention = paste(
      "suest2 repeats the full lnalpha cluster score on each observation;",
      "R uses one exact integrated-likelihood score per panel"
    )
  ),
  data = dat,
  balanced = list(
    model1 = balanced_model1,
    model2 = balanced_model2,
    prediction_mean = c(Y1 = 1.39122988319255, Y2 = 0.909426647884986),
    average_slope_x = c(Y1 = 0.506494047986875, Y2 = 0.221953416058067),
    system = balanced_system
  ),
  partial_higher = list(system = partial_system),
  disjoint = list(system = disjoint_system)
)

saveRDS(reference,
  "../../tests/testthat/fixtures/panel-poisson-re-stata.rds")
cat("PANEL_POISSON_RE_STATA_FIXTURE_COMPLETE=1\n")
