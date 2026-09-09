version 16
set more off
* Run after extracting all files into your working directory.
* No ado files are changed. Each benchmark writes its own text log.
local failures 0
foreach test in panel_fe panel_be_re panel_ml endogenous gee {
    capture noisily do suest_r_`test'_crosslang_benchmark.do
    local rc = _rc
    if `rc' {
        local failures = `failures' + 1
        di as error "BENCHMARK_FAILED: `test', return code `rc'"
    }
}
capture noisily do suest_r_pweight_lnvar_scale_diagnostic.do
local rc = _rc
if `rc' {
    local failures = `failures' + 1
    di as error "BENCHMARK_FAILED: pweight_lnvar_scale, return code `rc'"
}
di as result "BENCHMARK_BATCH_FAILURES=`failures'"
di as txt "Please upload all six benchmark/diagnostic .log files."
