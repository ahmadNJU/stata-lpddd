*==============================================================================================
* lpddd_example.do     example and tests for lpddd 0.1.0 (22 November 2026)
*
* Run in Stata:   do lpddd_example.do
* Needs:          lpddd, reghdfe, ftools  (lpdid is used in test 5 if it is installed)
*
* The data are simulated: 400 products (groups), each imported from China (eligible unit)
* and from all other suppliers (ineligible unit), 96 months. Three cohorts of products get
* a tariff on China from month 50, 53 or 56; 100 products never do. The true effect on
* imports from China is 0 in months 0 to 2 and -0.5 log points from month 3.
* Every product has its own level and trend, China has its own seasonal pattern and trend,
* and the taxed products have their own trend: a difference-in-differences is biased by
* these trends, the triple difference is not.
*
* Tests 2 and 3 stop with an error if lpddd does not reproduce the hand-made estimate.
*==============================================================================================
version 14
clear all
set more off
set seed 20261122

*----------------------------------------------------------------------------------------------
* Simulated data
*----------------------------------------------------------------------------------------------
set obs 400
generate int product = _n
generate byte chapter = 1 + mod(_n, 5)
generate int cohort = 0
replace cohort = 50 if _n <= 100
replace cohort = 53 if _n > 100 & _n <= 160
replace cohort = 56 if _n > 160 & _n <= 300
generate double a_g = rnormal(0, 1)
generate double trend_g = rnormal(0, 0.002)
generate double w = exp(a_g)
expand 2
bysort product: generate byte china = (_n == 1)
generate long unit = product * 10 + china
generate double a_u = rnormal(0, 0.5)
expand 96
bysort unit: generate int t = _n
generate byte month = 1 + mod(t - 1, 12)
generate int rel = t - cohort if cohort > 0
generate double tau = cond(cohort > 0 & china == 1 & rel >= 3 & rel < ., -0.5, 0)
generate double lny = 5 + a_g + a_u + trend_g * t                      ///
    + (0.3 + 0.2 * china) * sin(2 * _pi * month / 12)                  ///
    + 0.10 * china * (t / 96) + 0.05 * (cohort > 0) * (t / 96)         ///
    + tau + rnormal(0, 0.2)
generate double y = exp(lny)
label variable y "Imports"
tempfile sim
save `sim'

*----------------------------------------------------------------------------------------------
* 1. Main use: month by month and spans, seasonal base (expected: about 0 before month 3,
*    about -0.5 afterwards; spans 3 to 11 and 12 to 23 about -0.5)
*----------------------------------------------------------------------------------------------
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) pre(24) post(24) spans(-24 -13 -12 -1 0 2 3 11 12 23) graph
ereturn list

*----------------------------------------------------------------------------------------------
* 2. One stack, one month: lpddd against the triple difference of the four cell means
*----------------------------------------------------------------------------------------------
lpddd y if cohort == 50 | cohort == 0, unit(unit) time(t) group(product) cohort(cohort) ///
    eligible(china) log season(12) post(12) notable
matrix R = e(results)
scalar auto = .
forvalues i = 1/`=rowsof(R)' {
    if R[`i', 1] == 0 & R[`i', 2] == 12 scalar auto = R[`i', 4]
}
preserve
keep if cohort == 50 | cohort == 0
xtset unit t
generate double ly = ln(y)
generate double dy = ly - L24.ly
keep if t == 62
generate byte T = (cohort == 50)
regress dy i.T##i.china
scalar manual = _b[1.T#1.china]
restore
display as result "Test 2: lpddd " auto "   cell means " manual
assert reldif(auto, manual) < 1e-8

*----------------------------------------------------------------------------------------------
* 3. Three stacks, one span, group and sector effects: lpddd against a stacked regression
*    built by hand (the specification of a tariff-line trade study)
*----------------------------------------------------------------------------------------------
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(12 23) gfe strata(chapter)
matrix R = e(results)
scalar auto_b = R[1, 4]
scalar auto_se = R[1, 5]

local k 0
foreach g in 50 53 56 {
    local ++k
    use `sim', clear
    keep if cohort == `g' | cohort == 0
    generate int r = t - `g'
    generate byte T = (cohort == `g')
    bysort unit: egen double num = total(y * inrange(r, 12, 23))
    bysort unit: egen double den = total(y * inrange(r, -12, -1))
    bysort unit: keep if _n == 1
    generate double dy = ln(num) - ln(den)
    generate int stack = `k'
    tempfile s`k'
    save `s`k''
}
use `s1', clear
append using `s2'
append using `s3'
generate byte D = T * china
reghdfe dy D, absorb(i.stack#i.product i.stack#i.china#i.chapter) vce(cluster product)
scalar manual_b = _b[D]
scalar manual_se = _se[D]
display as result "Test 3: lpddd " auto_b " (" auto_se ")   by hand " manual_b " (" manual_se ")"
assert reldif(auto_b, manual_b) < 1e-8
assert reldif(auto_se, manual_se) < 1e-6
use `sim', clear

*----------------------------------------------------------------------------------------------
* 4. Weights and explicit aggregation of the stacks
*----------------------------------------------------------------------------------------------
lpddd y [aweight = w], unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(3 11 12 23) gfe strata(chapter)
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(3 11 12 23) aggregate(cohort)
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(3 11 12 23) aggregate(equal)

*----------------------------------------------------------------------------------------------
* 5. Difference-in-differences
*    a. other suppliers alone (expected: about 0.05 * 12/96 per year from the trend of the
*       taxed products, no effect of the tariff)
*    b. China alone without season, beside lpdid with never-treated controls if installed
*       (the two are expected to agree; this test only displays them)
*----------------------------------------------------------------------------------------------
lpddd y if china == 0, unit(unit) time(t) group(product) cohort(cohort) dd ///
    log season(12) noevent spans(-12 -1 3 11 12 23)

keep if china == 1
generate double ly = ln(y)
lpddd ly, unit(unit) time(t) group(product) cohort(cohort) dd pre(5) post(10)
capture which lpdid
if _rc == 0 {
    generate byte treat = (cohort > 0 & t >= cohort)
    lpdid ly, unit(unit) time(t) treat(treat) pre_window(5) post_window(10) nevertreated nograph
}
else display as text "lpdid is not installed; comparison skipped (ssc install lpdid)"

display _newline as result "lpddd_example.do finished: tests 2 and 3 passed."
