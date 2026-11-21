*! lpddd 0.1.0  22nov2026
*! Stacked local-projection triple differences (and differences-in-differences)
*! with a seasonal base period and estimates for spans of periods
*! Requires reghdfe (Correia), available from SSC

program define lpddd, eclass
    version 14.0
    syntax varname(numeric) [if] [in] [aweight pweight fweight] ,      ///
        Unit(varname) Time(varname numeric) GROup(varname)             ///
        COHort(varname numeric)                                        ///
        [ ELigible(varname numeric) DD PRE(integer 0) POST(integer 0)  ///
          SEASon(integer 1) SPANs(numlist integer) noEVent LOG IHS     ///
          FILLzero CONTrol(string) GFE STRata(varname) ABSorb(string)  ///
          AGGregate(string) CLuster(varname) Level(cilevel) GRAph      ///
          SAVing(string asis) noTABle ]

    local cmdline `"lpddd `0'"'
    local y `varlist'

    * ------------------------------------------------------------------ checks
    capture which reghdfe
    if _rc {
        display as error "lpddd requires reghdfe. Type: ssc install reghdfe  and  ssc install ftools"
        exit 199
    }
    if "`dd'" == "" & "`eligible'" == "" {
        display as error "specify eligible() for a triple difference, or dd for a difference-in-differences"
        exit 198
    }
    if "`dd'" != "" & "`eligible'" != "" {
        display as error "eligible() and dd may not be combined"
        exit 198
    }
    if "`dd'" != "" & "`gfe'" != "" {
        display as error "gfe may not be combined with dd (treatment does not vary within a group)"
        exit 198
    }
    if `season' < 1 {
        display as error "season() must be 1 or larger"
        exit 198
    }
    if `pre' < 0 | `post' < 0 {
        display as error "pre() and post() must not be negative"
        exit 198
    }
    if "`log'" != "" & "`ihs'" != "" {
        display as error "log and ihs may not be combined"
        exit 198
    }
    if "`control'" == "" local control "never"
    if !inlist("`control'", "never", "notyet") {
        display as error "control() must be never or notyet"
        exit 198
    }
    if "`aggregate'" == "" local aggregate "ols"
    if !inlist("`aggregate'", "ols", "cohort", "equal") {
        display as error "aggregate() must be ols, cohort or equal"
        exit 198
    }
    if "`event'" == "noevent" & "`spans'" == "" {
        display as error "nothing to estimate: noevent requires spans()"
        exit 198
    }
    local transf "none"
    if "`log'" != "" local transf "log"
    if "`ihs'" != "" local transf "ihs"

    * saving(filename [, replace])
    local savefile ""
    local savereplace ""
    if `"`saving'"' != "" {
        _lpddd_parsesave `saving'
        local savefile `"`r(file)'"'
        local savereplace "`r(replace)'"
    }

    * ------------------------------------------------------------------ spans
    * span j: periods a_j to b_j, compared with its reference periods
    *   after treatment (a_j >= 0): the same seasons of the last cycle before treatment
    *   before treatment (b_j < 0): the same periods one or more full cycles earlier
    local nsp 0
    if "`spans'" != "" {
        local nw : word count `spans'
        if mod(`nw', 2) != 0 {
            display as error "spans() needs pairs of numbers: first and last period of each span"
            exit 198
        }
        local nsp = `nw' / 2
        forvalues j = 1/`nsp' {
            local i1 = 2 * `j' - 1
            local i2 = 2 * `j'
            local a`j' : word `i1' of `spans'
            local b`j' : word `i2' of `spans'
            if `a`j'' > `b`j'' {
                display as error "spans(): the first period of a span must not exceed its last period"
                exit 198
            }
            if `a`j'' < 0 & `b`j'' >= 0 {
                display as error "spans(): a span may not contain periods before and after treatment"
                exit 198
            }
            local len`j' = `b`j'' - `a`j'' + 1
            if `a`j'' >= 0 {
                local rtype`j' "list"
                local ref`j' ""
                forvalues e = `a`j''/`b`j'' {
                    local r = `e' - `season' * (floor(`e' / `season') + 1)
                    local ref`j' : list ref`j' | r
                }
                local rlen`j' : word count `ref`j''
            }
            else {
                local rtype`j' "range"
                local m = ceil(`len`j'' / `season')
                local ra`j' = `a`j'' - `season' * `m'
                local rb`j' = `b`j'' - `season' * `m'
                local rlen`j' = `len`j''
            }
        }
    }

    * event window needed in every stack
    local lo = -max(`pre', `season')
    local hi = `post'
    forvalues j = 1/`nsp' {
        if `b`j'' > `hi' local hi = `b`j''
        if `a`j'' < `lo' local lo = `a`j''
        if "`rtype`j''" == "range" {
            if `ra`j'' < `lo' local lo = `ra`j''
        }
    }

    * ------------------------------------------------------------------ sample
    marksample touse, novarlist
    markout `touse' `unit' `time' `group', strok
    if "`eligible'" != "" markout `touse' `eligible'
    if "`strata'" != ""   markout `touse' `strata', strok
    if "`cluster'" != ""  markout `touse' `cluster', strok
    local wv ""
    if "`weight'" != "" {
        tempvar wv
        quietly generate double `wv' `exp'
        markout `touse' `wv'
    }

    * variables named in absorb()
    local absv ""
    if `"`absorb'"' != "" {
        local tmp : subinstr local absorb "#" " ", all
        foreach tk of local tmp {
            local v = regexr("`tk'", "^[a-z0-9]*\.", "")
            capture confirm variable `v'
            if !_rc local absv `absv' `v'
        }
    }

    preserve
    quietly keep if `touse'
    if _N == 0 {
        display as error "no observations"
        exit 2000
    }
    capture describe __lp_*, varlist
    if !_rc {
        display as error "variables with names starting with __lp_ are not allowed in the data"
        exit 110
    }
    local keepv `y' `unit' `time' `group' `cohort' `eligible' `strata' `cluster' `wv' `absv'
    local keepv : list uniq keepv
    keep `keepv'

    quietly egen long __lp_u = group(`unit')
    quietly egen long __lp_g = group(`group')

    capture assert `time' == int(`time')
    if _rc {
        display as error "time() must contain integers"
        exit 459
    }
    capture assert `cohort' == int(`cohort') if !missing(`cohort')
    if _rc {
        display as error "cohort() must contain integers (the period of first treatment), 0 or missing for never treated"
        exit 459
    }
    capture isid __lp_u `time'
    if _rc {
        display as error "unit() and time() do not uniquely identify the observations"
        exit 459
    }
    quietly generate byte __lp_never = missing(`cohort') | `cohort' == 0
    quietly bysort __lp_g (`cohort'): generate byte __lp_bad = (`cohort'[1] != `cohort'[_N])
    quietly count if __lp_bad
    if r(N) > 0 {
        display as error "cohort() must be constant within group()"
        exit 459
    }
    drop __lp_bad

    if "`eligible'" != "" {
        capture assert inlist(`eligible', 0, 1)
        if _rc {
            display as error "eligible() must be 0 or 1"
            exit 459
        }
        quietly bysort __lp_u (`eligible'): generate byte __lp_bad = (`eligible'[1] != `eligible'[_N])
        quietly count if __lp_bad
        if r(N) > 0 {
            display as error "eligible() must be constant within unit()"
            exit 459
        }
        drop __lp_bad
        quietly generate byte __lp_Q = `eligible'
    }
    else quietly generate byte __lp_Q = 1

    local strk ""
    if "`strata'" != "" {
        quietly egen long __lp_str = group(`strata')
        local strk "__lp_str"
    }
    if "`cluster'" == "" local clv "__lp_g"
    else {
        quietly egen long __lp_cl = group(`cluster')
        local clv "__lp_cl"
    }
    local wkeep ""
    if "`weight'" != "" {
        quietly generate double __lp_w = `wv'
        quietly bysort __lp_u (__lp_w): generate byte __lp_bad = (__lp_w[1] != __lp_w[_N])
        quietly count if __lp_bad
        if r(N) > 0 {
            display as error "weights must be constant within unit()"
            exit 459
        }
        drop __lp_bad
        local wkeep "__lp_w"
    }

    quietly levelsof `cohort' if !__lp_never, local(glist)
    local ng : word count `glist'
    if `ng' == 0 {
        display as error "no treated cohort: cohort() is 0 or missing in all observations"
        exit 459
    }
    quietly count if __lp_never
    if r(N) == 0 & "`control'" == "never" {
        display as error "no never-treated group; consider control(notyet)"
        exit 459
    }

    * fixed effects
    if "`dd'" != "" {
        local fe "i.__lp_stack"
        if "`strata'" != "" local fe "i.__lp_stack#i.__lp_str"
    }
    else {
        if "`gfe'" != "" local fe1 "i.__lp_stack#i.__lp_g"
        else if "`strata'" != "" local fe1 "i.__lp_stack#i.__lp_T#i.__lp_str"
        else local fe1 "i.__lp_stack#i.__lp_T"
        if "`strata'" != "" local fe2 "i.__lp_stack#i.__lp_Q#i.__lp_str"
        else local fe2 "i.__lp_stack#i.__lp_Q"
        local fe "`fe1' `fe2'"
    }
    if `"`absorb'"' != "" local fe `"`fe' `absorb'"'

    local carry __lp_u __lp_g __lp_T __lp_Q __lp_D __lp_stack `clv' `wkeep' `strk' `absv'
    local carry : list uniq carry

    tempfile base
    quietly save `base'

    * ------------------------------------------------------------------ stacks
    local k 0
    foreach g of local glist {
        local ++k
        tempfile wk`k' ev`k' sp`k'
        quietly {
            use `base', clear
            if "`control'" == "never" keep if `cohort' == `g' | __lp_never
            else keep if `cohort' == `g' | __lp_never | (`cohort' > `g' + `hi' & !missing(`cohort'))
            generate long __lp_rel = `time' - `g'
            keep if inrange(__lp_rel, `lo', `hi')
            generate byte __lp_T = (`cohort' == `g')
            generate int __lp_stack = `k'
            if "`dd'" != "" generate byte __lp_D = __lp_T
            else generate byte __lp_D = __lp_T * __lp_Q
            count if __lp_D == 1
            local ntr`k' = r(N)
            save `wk`k''

            * spans: one observation per unit
            if `nsp' > 0 {
                forvalues j = 1/`nsp' {
                    generate byte __lp_i1 = inrange(__lp_rel, `a`j'', `b`j'')
                    if "`rtype`j''" == "range" generate byte __lp_i0 = inrange(__lp_rel, `ra`j'', `rb`j'')
                    else {
                        generate byte __lp_i0 = 0
                        foreach r of local ref`j' {
                            replace __lp_i0 = 1 if __lp_rel == `r'
                        }
                    }
                    bysort __lp_u: egen double __lp_s1 = total(`y' * __lp_i1)
                    bysort __lp_u: egen double __lp_s0 = total(`y' * __lp_i0)
                    bysort __lp_u: egen long __lp_n1 = total(__lp_i1 * !missing(`y'))
                    bysort __lp_u: egen long __lp_n0 = total(__lp_i0 * !missing(`y'))
                    replace __lp_s1 = __lp_s1 / `len`j''
                    replace __lp_s0 = __lp_s0 / `rlen`j''
                    if "`transf'" == "log" generate double __lp_sp`j' = ln(__lp_s1) - ln(__lp_s0) if __lp_s1 > 0 & __lp_s0 > 0
                    else if "`transf'" == "ihs" generate double __lp_sp`j' = asinh(__lp_s1) - asinh(__lp_s0)
                    else generate double __lp_sp`j' = __lp_s1 - __lp_s0
                    if "`fillzero'" == "" replace __lp_sp`j' = . if __lp_n1 < `len`j'' | __lp_n0 < `rlen`j''
                    drop __lp_i1 __lp_i0 __lp_s1 __lp_s0 __lp_n1 __lp_n0
                }
                bysort __lp_u (`time'): keep if _n == 1
                keep `carry' __lp_sp*
                save `sp`k''
                use `wk`k'', clear
            }

            * event time: one observation per unit and period
            if "`event'" != "noevent" {
                if "`transf'" == "log" generate double __lp_f = ln(`y') if `y' > 0
                else if "`transf'" == "ihs" generate double __lp_f = asinh(`y')
                else generate double __lp_f = `y'
                xtset __lp_u `time'
                generate long __lp_d = `season' * (floor(__lp_rel / `season') + 1)
                generate double __lp_dy = .
                levelsof __lp_d if inrange(__lp_rel, -`pre', `post') & __lp_d != 0, local(dl)
                foreach d of local dl {
                    if `d' > 0 {
                        replace __lp_dy = __lp_f - L`d'.__lp_f if __lp_d == `d'
                    }
                    else {
                        local ad = -`d'
                        replace __lp_dy = __lp_f - F`ad'.__lp_f if __lp_d == `d'
                    }
                }
                keep if inrange(__lp_rel, -`pre', `post') & __lp_d != 0 & !missing(__lp_dy)
                keep `carry' __lp_rel __lp_dy
                save `ev`k''
            }
        }
    }

    * ------------------------------------------------------------------ estimates
    * columns of the results matrix:
    *   type (0 event time, 1 span), from, to, b, se, t, p, ll, ul, N, ref (1 = reference period)
    tempname R bb vv sb sse sN sdf st sp sll sul crit
    local names ""
    local Nmax 0
    local wopt ""
    if "`weight'" != "" local wopt "wtype(`weight')"

    if "`event'" != "noevent" {
        quietly {
            use `ev1', clear
            forvalues k = 2/`ng' {
                append using `ev`k''
            }
            if "`aggregate'" != "ols" {
                forvalues k = 1/`ng' {
                    generate byte __lp_D`k' = __lp_D * (__lp_stack == `k')
                }
            }
        }
        forvalues e = -`pre'/`post' {
            local d = `season' * (floor(`e' / `season') + 1)
            if `d' == 0 {
                matrix `R' = nullmat(`R') \ (0, `e', `e', 0, ., ., ., 0, 0, ., 1)
                continue
            }
            _lpddd_fit __lp_dy if __lp_rel == `e', d(__lp_D) fe(`fe') clv(`clv') ///
                agg(`aggregate') nst(`ng') `wopt'
            scalar `sb'  = r(b)
            scalar `sse' = r(se)
            scalar `sN'  = r(N)
            scalar `sdf' = r(df)
            scalar `st'  = .
            scalar `sp'  = .
            scalar `sll' = .
            scalar `sul' = .
            if `sse' > 0 & `sse' < . {
                scalar `st'   = `sb' / `sse'
                scalar `sp'   = 2 * ttail(`sdf', abs(`st'))
                scalar `crit' = invttail(`sdf', (100 - `level') / 200)
                scalar `sll'  = `sb' - `crit' * `sse'
                scalar `sul'  = `sb' + `crit' * `sse'
                if `e' < 0 {
                    local ae = -`e'
                    local nm "m`ae'"
                }
                else local nm "p`e'"
                matrix `bb' = nullmat(`bb'), `sb'
                matrix `vv' = nullmat(`vv'), `sse'^2
                local names `names' `nm'
                if `sN' > `Nmax' local Nmax = `sN'
            }
            matrix `R' = nullmat(`R') \ (0, `e', `e', `sb', `sse', `st', `sp', `sll', `sul', `sN', 0)
        }
    }

    if `nsp' > 0 {
        quietly {
            use `sp1', clear
            forvalues k = 2/`ng' {
                append using `sp`k''
            }
            if "`aggregate'" != "ols" {
                forvalues k = 1/`ng' {
                    generate byte __lp_D`k' = __lp_D * (__lp_stack == `k')
                }
            }
        }
        forvalues j = 1/`nsp' {
            _lpddd_fit __lp_sp`j', d(__lp_D) fe(`fe') clv(`clv') ///
                agg(`aggregate') nst(`ng') `wopt'
            scalar `sb'  = r(b)
            scalar `sse' = r(se)
            scalar `sN'  = r(N)
            scalar `sdf' = r(df)
            scalar `st'  = .
            scalar `sp'  = .
            scalar `sll' = .
            scalar `sul' = .
            if `sse' > 0 & `sse' < . {
                scalar `st'   = `sb' / `sse'
                scalar `sp'   = 2 * ttail(`sdf', abs(`st'))
                scalar `crit' = invttail(`sdf', (100 - `level') / 200)
                scalar `sll'  = `sb' - `crit' * `sse'
                scalar `sul'  = `sb' + `crit' * `sse'
                matrix `bb' = nullmat(`bb'), `sb'
                matrix `vv' = nullmat(`vv'), `sse'^2
                local names `names' span`j'
                if `sN' > `Nmax' local Nmax = `sN'
            }
            matrix `R' = nullmat(`R') \ (1, `a`j'', `b`j'', `sb', `sse', `st', `sp', `sll', `sul', `sN', 0)
        }
    }
    matrix colnames `R' = type from to b se t p ll ul N ref

    * ------------------------------------------------------------------ results as data: graph, saving
    if "`graph'" != "" | `"`savefile'"' != "" {
        quietly {
            drop _all
            svmat double `R', names(col)
            label variable type "0 event time, 1 span"
            label variable from "First period"
            label variable to   "Last period"
            label variable b    "Estimate"
            label variable se   "Standard error"
            label variable ll   "Lower limit"
            label variable ul   "Upper limit"
            label variable N    "Observations"
            label variable ref  "1 = reference period (zero by construction)"
        }
        if `"`savefile'"' != "" {
            quietly save `savefile', `savereplace'
        }
        if "`graph'" != "" {
            quietly count if type == 0
            if r(N) > 0 {
                twoway (rarea ll ul from if type == 0, color(gs13) lwidth(none))        ///
                       (line b from if type == 0, lcolor(black) lwidth(medthick)),      ///
                       yline(0, lcolor(gs8) lwidth(thin))                               ///
                       xline(-0.5, lcolor(gs8) lpattern(shortdash))                     ///
                       xtitle("Periods from the first period of treatment")             ///
                       ytitle("Estimate")                                               ///
                       legend(order(2 "Estimate" 1 "`level' percent interval") rows(1)) ///
                       graphregion(color(white))
            }
        }
    }
    restore

    * ------------------------------------------------------------------ display
    local design "triple difference (treated group x eligible unit)"
    if "`dd'" != "" local design "difference-in-differences (treated group)"
    display _newline as text "Stacked local projections: `design'"
    display as text "Outcome: " as result "`y'" as text "   transformation: " as result "`transf'" ///
        as text "   seasonal cycle: " as result "`season'"
    display as text "Stacks (treated cohorts): " as result "`ng'" as text "   controls: " as result "`control'" ///
        as text "   aggregation: " as result "`aggregate'"
    local k 0
    foreach g of local glist {
        local ++k
        display as text "  stack `k': first treated period " as result "`g'" ///
            as text ", treated observations in the window " as result "`ntr`k''"
    }
    display as text "Fixed effects: " as result `"`fe'"'
    display as text "Standard errors clustered by: " as result cond("`cluster'" == "", "`group'", "`cluster'")

    if "`table'" != "notable" {
        local nr = rowsof(`R')
        display _newline as text "{hline 15}{c TT}{hline 66}"
        display as text %14s "Periods" " {c |}" %11s "Estimate" %11s "Std. err." %9s "t" %9s "P>|t|" ///
            %12s "[`level'% conf." %11s "interval]" %9s "Obs."
        display as text "{hline 15}{c +}{hline 66}"
        forvalues i = 1/`nr' {
            if `R'[`i', 1] == 0 local lab = string(`R'[`i', 2])
            else local lab = string(`R'[`i', 2]) + " to " + string(`R'[`i', 3])
            if `R'[`i', 11] == 1 {
                display as text %14s "`lab'" " {c |}" as text %11s "0" as text "   (reference)"
            }
            else {
                display as text %14s "`lab'" " {c |}" as result %11.4f `R'[`i', 4] %11.4f `R'[`i', 5] ///
                    %9.2f `R'[`i', 6] %9.3f `R'[`i', 7] %12.4f `R'[`i', 8] %11.4f `R'[`i', 9] %9.0f `R'[`i', 10]
            }
            if `i' < `nr' {
                if `R'[`i', 1] == 0 & `R'[`i' + 1, 1] == 1 display as text "{hline 15}{c +}{hline 66}"
            }
        }
        display as text "{hline 15}{c BT}{hline 66}"
    }

    * ------------------------------------------------------------------ stored results
    ereturn clear
    if "`names'" != "" {
        tempname VV
        matrix `VV' = diag(`vv')
        matrix colnames `bb' = `names'
        matrix colnames `VV' = `names'
        matrix rownames `VV' = `names'
        ereturn post `bb' `VV', obs(`Nmax') depname(`y')
    }
    ereturn matrix results = `R'
    ereturn scalar n_stacks = `ng'
    ereturn scalar season = `season'
    ereturn scalar pre = `pre'
    ereturn scalar post = `post'
    ereturn scalar level = `level'
    ereturn local cohorts "`glist'"
    ereturn local absorb `"`fe'"'
    ereturn local aggregate "`aggregate'"
    ereturn local control "`control'"
    ereturn local transform "`transf'"
    ereturn local design = cond("`dd'" != "", "dd", "ddd")
    ereturn local clustvar = cond("`cluster'" == "", "`group'", "`cluster'")
    ereturn local depvar "`y'"
    ereturn local cmdline `"`cmdline'"'
    ereturn local cmd "lpddd"
end


* one regression: the coefficient on the treatment indicator, or an explicitly
* weighted average of stack-specific coefficients estimated in one regression
program define _lpddd_fit, rclass
    version 14.0
    syntax varname [if] , d(varname) fe(string) clv(varname) agg(string) nst(integer) [ wtype(string) ]
    local wx ""
    if "`wtype'" != "" local wx "[`wtype' = __lp_w]"
    tempname b se N df
    scalar `b' = .
    scalar `se' = .
    scalar `N' = .
    scalar `df' = .
    if "`agg'" == "ols" {
        capture reghdfe `varlist' `d' `if' `wx', absorb(`fe') vce(cluster `clv')
        if _rc == 0 {
            scalar `N' = e(N)
            scalar `df' = e(df_r)
            capture scalar `b' = _b[`d']
            if _rc == 0 scalar `se' = _se[`d']
        }
    }
    else {
        local dl ""
        forvalues k = 1/`nst' {
            local dl `dl' __lp_D`k'
        }
        capture reghdfe `varlist' `dl' `if' `wx', absorb(`fe') vce(cluster `clv')
        if _rc == 0 {
            scalar `N' = e(N)
            scalar `df' = e(df_r)
            tempvar es
            quietly generate byte `es' = e(sample)
            local tot 0
            forvalues k = 1/`nst' {
                local w`k' 0
                tempname sek
                capture scalar `sek' = _se[__lp_D`k']
                if _rc == 0 {
                    if `sek' > 0 & `sek' < . {
                        if "`agg'" == "equal" local w`k' 1
                        else {
                            quietly count if `es' & __lp_D`k' == 1
                            local w`k' = r(N)
                        }
                    }
                }
                local tot = `tot' + `w`k''
            }
            if `tot' > 0 {
                local expr ""
                forvalues k = 1/`nst' {
                    if `w`k'' > 0 {
                        local c = `w`k'' / `tot'
                        if "`expr'" == "" local expr "`c' * _b[__lp_D`k']"
                        else local expr "`expr' + `c' * _b[__lp_D`k']"
                    }
                }
                capture lincom `expr'
                if _rc == 0 {
                    scalar `b' = r(estimate)
                    scalar `se' = r(se)
                }
            }
        }
    }
    return scalar b = `b'
    return scalar se = `se'
    return scalar N = `N'
    return scalar df = `df'
end


program define _lpddd_parsesave, rclass
    version 14.0
    syntax anything(name=fn) [, replace]
    return local file `"`fn'"'
    return local replace "`replace'"
end
