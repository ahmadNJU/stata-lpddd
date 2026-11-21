*==============================================================================================
* lpddd_example.do     example and tests for lpddd 0.4.0 (22 November 2026)
*
* Run in Stata:   do lpddd_example.do
* Needs:          lpddd. Tests 3 and 6 need reghdfe and ftools and are skipped without them;
*                 test 5 uses lpdid if it is installed. Run with variable abbreviation off.
*
* The data are simulated: 400 products (groups), each imported from China (eligible unit)
* and from all other suppliers (ineligible unit), 96 months. Three cohorts of products get
* a tariff on China from month 50, 53 or 56; 100 products never do. The true effect on
* imports from China is 0 in months 0 to 2 and -0.5 log points from month 3.
* Every product has its own level and trend, China has its own seasonal pattern and trend,
* and the taxed products have their own trend: a difference-in-differences is biased by
* these trends, the triple difference is not.
*
* Tests 2, 3, 6, 7, 8, 9 and 10 stop with an error if lpddd does not reproduce the estimate it is
* compared with.
*==============================================================================================
version 15.1
clear all
set more off
set varabbrev off
capture which reghdfe
local has_reghdfe = (_rc == 0)

* To test lpddd without reghdfe, change the 0 in the next line to 1: every command then
* uses the built-in commands (regress and xtreg) and tests 3 and 6 are skipped.
local builtin_only 0
local eng ""
if `builtin_only' {
    local has_reghdfe 0
    local eng "engine(builtin)"
}
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
    log season(12) pre(24) post(24) spans(-24 -13 -12 -1 0 2 3 11 12 23) graph `eng'
ereturn list

*----------------------------------------------------------------------------------------------
* 2. One stack, one month: lpddd against the triple difference of the four cell means
*----------------------------------------------------------------------------------------------
lpddd y if cohort == 50 | cohort == 0, unit(unit) time(t) group(product) cohort(cohort) ///
    eligible(china) log season(12) post(12) notable `eng'
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
if `has_reghdfe' {
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(12 23) gfe strata(chapter) engine(reghdfe)
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
}
else display as text "reghdfe is not installed; test 3 skipped"
use `sim', clear

*----------------------------------------------------------------------------------------------
* 4. Weights and explicit aggregation of the stacks
*----------------------------------------------------------------------------------------------
lpddd y [aweight = w], unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(3 11 12 23) gfe strata(chapter) `eng'
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(3 11 12 23) aggregate(cohort) `eng'
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(3 11 12 23) aggregate(equal) `eng'

*----------------------------------------------------------------------------------------------
* 5. Difference-in-differences
*    a. other suppliers alone (expected: about 0.05 * 12/96 per year from the trend of the
*       taxed products, no effect of the tariff)
*    b. China alone without season, beside lpdid with never-treated controls if installed
*       (the two are expected to agree; this test only displays them)
*----------------------------------------------------------------------------------------------
lpddd y if china == 0, unit(unit) time(t) group(product) cohort(cohort) dd ///
    log season(12) noevent spans(-12 -1 3 11 12 23) `eng'

keep if china == 1
generate double ly = ln(y)
lpddd ly, unit(unit) time(t) group(product) cohort(cohort) dd pre(5) post(10) `eng'
capture which lpdid
if _rc == 0 {
    generate byte treat = (cohort > 0 & t >= cohort)
    lpdid ly, unit(unit) time(t) treat(treat) pre_window(5) post_window(10) nevertreated nograph
}
else display as text "lpdid is not installed; comparison skipped (ssc install lpdid)"

use `sim', clear

*----------------------------------------------------------------------------------------------
* 6. Built-in commands against reghdfe for five specifications: the coefficients and the
*    standard errors must agree
*----------------------------------------------------------------------------------------------
if `has_reghdfe' {
    forvalues s = 1/5 {
        if `s' == 1 local spec ""
        if `s' == 2 local spec "gfe strata(chapter)"
        if `s' == 3 local spec "strata(chapter) aggregate(cohort)"
        if `s' == 4 local spec "gfe"
        if `s' == 5 local spec "gfe strata(chapter) aggregate(cohort)"
        lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
            log season(12) noevent spans(12 23) `spec' engine(reghdfe) notable
        matrix A = e(results)
        lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
            log season(12) noevent spans(12 23) `spec' engine(builtin) notable
        matrix B = e(results)
        display as result "Test 6, specification `s' (`spec'): reghdfe " A[1, 4] " (" A[1, 5] ")   built-in " B[1, 4] " (" B[1, 5] ")"
        display as result "        relative difference: coefficient " reldif(A[1, 4], B[1, 4]) ", standard error " abs(A[1, 5] - B[1, 5]) / abs(A[1, 5])
        assert reldif(A[1, 4], B[1, 4]) < 1e-8
        assert abs(A[1, 5] - B[1, 5]) / abs(A[1, 5]) < 1e-6
    }
}
else display as text "reghdfe is not installed; test 6 skipped"

*----------------------------------------------------------------------------------------------
* 7. Explicit aggregation: with equal weights the estimate must be the simple mean of the
*    three estimates obtained stack by stack, and e(stack_b) must hold these three estimates
*    (built-in commands, which solve the regression exactly)
*----------------------------------------------------------------------------------------------
scalar avg = 0
local k 0
foreach g in 50 53 56 {
    local ++k
    lpddd y if cohort == `g' | cohort == 0, unit(unit) time(t) group(product) cohort(cohort) ///
        eligible(china) log season(12) noevent spans(12 23) engine(builtin) notable
    matrix A = e(results)
    scalar one`k' = A[1, 4]
    scalar avg = avg + A[1, 4] / 3
}
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(12 23) aggregate(equal) engine(builtin)
matrix B = e(results)
matrix S = e(stack_b)
matrix W = e(stack_w)
matrix list S
matrix list W
display as result "Test 7: mean of the three stacks " %18.15f avg "   aggregate(equal) " %18.15f B[1, 4]
assert reldif(avg, B[1, 4]) < 1e-10
forvalues k = 1/3 {
    assert reldif(one`k', S[1, `k']) < 1e-10
}
ereturn list

*----------------------------------------------------------------------------------------------
* 8. Covariances between periods, joint test and e(sample)
*    a. the diagonal of e(V) must hold the squared standard errors of the table
*    b. the correlation between the estimates for months 3 and 4 must equal the one from a
*       single regression on both months, with every term interacted with the month
*    c. test and lincom must run after lpddd; the joint test of the months before
*       treatment is displayed
*    d. with nocov, e(V) is diagonal and the estimates are unchanged
*----------------------------------------------------------------------------------------------
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) pre(24) post(6) spans(3 6) `eng'
matrix V = e(V)
matrix R = e(results)
matrix bfull = e(b)
display as result "Test 8: e(V) is " e(vfull) ";  joint test before treatment F(" e(df_pre) ", " e(df_r) ") = " e(F_pre) ", p = " e(p_pre)
assert "`e(vfull)'" == "full"
forvalues i = 1/`=rowsof(R)' {
    if R[`i', 1] == 0 & R[`i', 2] == 3 scalar se3 = R[`i', 5]
    if R[`i', 1] == 0 & R[`i', 2] == 4 scalar se4 = R[`i', 5]
}
local j3 = colnumb(V, "p3")
local j4 = colnumb(V, "p4")
assert reldif(V[`j3', `j3'], se3^2) < 1e-10
assert reldif(V[`j4', `j4'], se4^2) < 1e-10
scalar corr_auto = V[`j3', `j4'] / sqrt(V[`j3', `j3'] * V[`j4', `j4'])
quietly count if e(sample)
display as result "        observations marked by e(sample): " r(N) " of " _N
assert r(N) > 0
test p3 = p4
lincom (p3 + p4 + p5) / 3
test m24 m23 m22 m21 m20 m19 m18 m17 m16 m15 m14 m13
display as result "        test command, months before treatment: F = " r(F) ", p = " r(p)
assert reldif(r(F), e(F_pre)) < 1e-8

local k 0
foreach g in 50 53 56 {
    local ++k
    use `sim', clear
    keep if cohort == `g' | cohort == 0
    xtset unit t
    generate double ly = ln(y)
    generate double dy = ly - L12.ly
    generate int h = t - `g'
    keep if h == 3 | h == 4
    generate byte T = (cohort == `g')
    generate int stack = `k'
    tempfile c`k'
    save `c`k''
}
use `c1', clear
append using `c2'
append using `c3'
generate byte D = T * china
generate byte D3 = D * (h == 3)
generate byte D4 = D * (h == 4)
regress dy D3 D4 i.h#i.stack#i.T i.h#i.stack#i.china, vce(cluster product)
matrix VM = e(V)
scalar corr_manual = VM[1, 2] / sqrt(VM[1, 1] * VM[2, 2])
display as result "Test 8: correlation of the estimates for months 3 and 4: lpddd " corr_auto "   one regression " corr_manual
assert reldif(corr_auto, corr_manual) < 1e-8
use `sim', clear

lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) pre(24) post(6) spans(3 6) nocov notable `eng'
matrix V0 = e(V)
matrix b0 = e(b)
assert "`e(vfull)'" == "diagonal"
assert V0[`j3', `j4'] == 0
assert reldif(V0[`j3', `j3'], se3^2) < 1e-10
assert mreldif(b0, bfull) < 1e-12

*----------------------------------------------------------------------------------------------
* 9. State dependence with by(): the state is the group of chapters (0 = chapters 1 and 2,
*    1 = chapters 3 to 5); the true effect is the same in both states
*    a. the estimates by state must equal those of the command run on each state alone
*    b. the test of equal effects must equal the one of the test command
*    c. built-in commands and reghdfe must agree
*----------------------------------------------------------------------------------------------
generate byte big = (chapter >= 3)
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) pre(24) post(6) spans(12 23) by(big) graph `eng'
matrix RB = e(results)
matrix H = e(het)
matrix list H
ereturn list
forvalues i = 1/`=rowsof(RB)' {
    if RB[`i', 1] == 1 & RB[`i', 12] == 1 scalar by1 = RB[`i', 4]
    if RB[`i', 1] == 1 & RB[`i', 12] == 2 scalar by2 = RB[`i', 4]
}
* the observations of the two states add up to those of the regression
scalar nstate = 0
forvalues i = 1/`=rowsof(RB)' {
    if RB[`i', 1] == 1 scalar nstate = nstate + RB[`i', 10]
}
display as result "Test 9: observations of the two states together " nstate ",  of the regression " e(N)
assert nstate == e(N)
test span1_1 = span1_2
scalar F_test = r(F)
assert reldif(F_test, H[rowsof(H), 4]) < 1e-8
lpddd y if big == 0, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(12 23) notable `eng'
matrix A = e(results)
scalar alone1 = A[1, 4]
lpddd y if big == 1, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) noevent spans(12 23) notable `eng'
matrix A = e(results)
scalar alone2 = A[1, 4]
display as result "Test 9: state 1, by() " by1 "   alone " alone1 ";   state 2, by() " by2 "   alone " alone2
assert reldif(by1, alone1) < 1e-8
assert reldif(by2, alone2) < 1e-8
if `has_reghdfe' {
    lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
        log season(12) noevent spans(3 11 12 23) gfe strata(chapter) by(big) engine(reghdfe) notable
    matrix A = e(results)
    matrix HA = e(het)
    lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
        log season(12) noevent spans(3 11 12 23) gfe strata(chapter) by(big) engine(builtin) notable
    matrix B = e(results)
    matrix HB = e(het)
    forvalues i = 1/`=rowsof(A)' {
        display as result "Test 9, row `i': reghdfe " A[`i', 4] " (" A[`i', 5] ")   built-in " B[`i', 4] " (" B[`i', 5] ")"
        assert reldif(A[`i', 4], B[`i', 4]) < 1e-8
        assert abs(A[`i', 5] - B[`i', 5]) / abs(A[`i', 5]) < 1e-6
    }
    assert mreldif(HA, HB) < 1e-6
}

*----------------------------------------------------------------------------------------------
* 10. Coding of the dates and other frequencies
*    a. months coded as yyyymm: with tcode(yyyymm) the estimates must equal those with the
*       count of months; without tcode() the command must stop with error 451
*    b. quarterly data (the months added up to quarters), season(4): against the triple
*       difference of the four cell means; quarters coded as yyyyq with tcode(yyyyq)
*    c. annual data (the months added up to years), no season: against the cell means
*----------------------------------------------------------------------------------------------
use `sim', clear
generate long tym = (2011 + floor((t - 1) / 12)) * 100 + month
generate long cym = 0
replace cym = (2011 + floor((cohort - 1) / 12)) * 100 + 1 + mod(cohort - 1, 12) if cohort > 0
lpddd y, unit(unit) time(t) group(product) cohort(cohort) eligible(china) ///
    log season(12) pre(24) post(6) spans(3 6) notable `eng'
matrix b_count = e(b)
lpddd y, unit(unit) time(tym) group(product) cohort(cym) eligible(china) ///
    log season(12) pre(24) post(6) spans(3 6) tcode(yyyymm) notable `eng'
matrix b_code = e(b)
quietly count if e(sample)
scalar n_es = r(N)
display as result "Test 10a: largest difference between tcode(yyyymm) and the count of months " mreldif(b_count, b_code) ";  e(sample) " n_es
assert mreldif(b_count, b_code) < 1e-12
assert n_es == 26000
capture noisily lpddd y, unit(unit) time(tym) group(product) cohort(cym) eligible(china) ///
    log season(12) pre(24) post(6) notable `eng'
display as result "Test 10a: return code without tcode() " _rc " (451 expected)"
assert _rc == 451

* quarterly
use `sim', clear
generate int q = ceil(t / 3)
generate int cq = cond(cohort > 0, ceil(cohort / 3), 0)
collapse (sum) y, by(unit product china chapter q cq)
lpddd y, unit(unit) time(q) group(product) cohort(cq) eligible(china) ///
    log season(4) pre(8) post(8) spans(4 7) `eng'
matrix b_q = e(b)
lpddd y if cq == 17 | cq == 0, unit(unit) time(q) group(product) cohort(cq) eligible(china) ///
    log season(4) post(4) notable `eng'
matrix R = e(results)
scalar auto_q = .
forvalues i = 1/`=rowsof(R)' {
    if R[`i', 1] == 0 & R[`i', 2] == 1 scalar auto_q = R[`i', 4]
}
preserve
keep if cq == 17 | cq == 0
xtset unit q
generate double ly = ln(y)
generate double dy = ly - L4.ly
keep if q == 18
generate byte T = (cq == 17)
quietly regress dy i.T##i.china
scalar manual_q = _b[1.T#1.china]
restore
display as result "Test 10b: quarter 1 after treatment, lpddd " auto_q "   cell means " manual_q
assert reldif(auto_q, manual_q) < 1e-8
generate long qcode = (2011 + floor((q - 1) / 4)) * 10 + 1 + mod(q - 1, 4)
generate long cqcode = cond(cq > 0, (2011 + floor((cq - 1) / 4)) * 10 + 1 + mod(cq - 1, 4), 0)
lpddd y, unit(unit) time(qcode) group(product) cohort(cqcode) eligible(china) ///
    log season(4) pre(8) post(8) spans(4 7) tcode(yyyyq) notable `eng'
matrix b_qcode = e(b)
assert mreldif(b_q, b_qcode) < 1e-12
capture noisily lpddd y, unit(unit) time(qcode) group(product) cohort(cqcode) eligible(china) ///
    log season(4) pre(8) post(8) notable `eng'
assert _rc == 451

* annual
use `sim', clear
generate int yr = 2010 + ceil(t / 12)
generate int cy = cond(cohort > 0, 2010 + ceil(cohort / 12), 0)
collapse (sum) y, by(unit product china chapter yr cy)
lpddd y, unit(unit) time(yr) group(product) cohort(cy) eligible(china) ///
    log pre(3) post(3) spans(1 3) `eng'
matrix R = e(results)
scalar auto_y = .
forvalues i = 1/`=rowsof(R)' {
    if R[`i', 1] == 0 & R[`i', 2] == 1 scalar auto_y = R[`i', 4]
}
preserve
xtset unit yr
generate double ly = ln(y)
generate double dy = ly - L2.ly
keep if yr == 2016
generate byte T = (cy == 2015)
quietly regress dy i.T##i.china
scalar manual_y = _b[1.T#1.china]
restore
display as result "Test 10c: year 1 after treatment, lpddd " auto_y "   cell means " manual_y
assert reldif(auto_y, manual_y) < 1e-8

display _newline as result "lpddd_example.do finished: all tests that were run passed."
