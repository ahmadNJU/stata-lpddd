*! lpddd 0.1.1  22nov2026
*! Stacked local-projection triple differences (and differences-in-differences)
*! with a seasonal base period and estimates for spans of periods
*! Uses reghdfe (Correia, SSC) if it is installed; otherwise regress and xtreg

program define lpddd, eclass
    version 15.1
    syntax varname(numeric) [if] [in] [aweight pweight fweight] ,      ///
        Unit(varname) Time(varname numeric) GROup(varname)             ///
        COHort(varname numeric)                                        ///
        [ ELigible(varname numeric) DD PRE(integer 0) POST(integer 0)  ///
          SEASon(integer 1) SPANs(numlist integer) noEVent LOG IHS     ///
          FILLzero CONTrol(string) GFE STRata(varname) ABSorb(string)  ///
          AGGregate(string) CLuster(varname) Level(cilevel) GRAph      ///
          ENGine(string)                                               ///
          SAVing(string asis) noTABle ]

    local cmdline `"lpddd `0'"'
    local y `varlist'

    * ------------------------------------------------------------------ checks
    if !inlist("`engine'", "", "reghdfe", "builtin") {
        display as error "engine() must be reghdfe or builtin"
        exit 198
    }
    capture which reghdfe
    local hasreghdfe = (_rc == 0)
    if "`engine'" == "" {
        if `hasreghdfe' local engine "reghdfe"
        else local engine "builtin"
    }
    if "`engine'" == "reghdfe" & !`hasreghdfe' {
        display as error "engine(reghdfe) requires reghdfe. Type: ssc install reghdfe  and  ssc install ftools"
        exit 199
    }
    if "`engine'" == "builtin" & `"`absorb'"' != "" {
        display as error "absorb() requires reghdfe (ssc install reghdfe; ssc install ftools)"
        exit 198
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
        exit 451
    }
    capture assert `cohort' == int(`cohort') if !missing(`cohort')
    if _rc {
        display as error "cohort() must contain integers (the period of first treatment), 0 or missing for never treated"
        exit 459
    }
    capture isid __lp_u `time'
    if _rc {
        display as error "repeated time values within unit(): unit() and time() must identify the observations"
        exit 451
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
            exit 450
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
        exit 2001
    }
    quietly count if __lp_never
    if r(N) == 0 & "`control'" == "never" {
        display as error "no never-treated group; consider control(notyet)"
        exit 2001
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
    * built-in commands: "one" = one set of effects (regress), "two" = two sets (regress),
    * "xt" = group effects absorbed by xtreg, fe and the other set as indicators
    local mode "two"
    if "`dd'" != "" local mode "one"
    else if "`gfe'" != "" local mode "xt"

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
    tempname R bb vv sb sse sN sdf sncl st sp sll sul crit BS WS tmpm
    local names ""
    local Nmax 0
    local Nclmax 0
    local lastrc 0
    local nfail 0
    local fitopt d(__lp_D) fe(`fe') clv(`clv') agg(`aggregate') nst(`ng') engine(`engine') mode(`mode')
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
            if "`engine'" == "builtin" {
                if "`mode'" == "one" egen long __lp_fa = group(__lp_stack `strk')
                else if "`mode'" == "xt" egen long __lp_fa = group(__lp_stack __lp_g)
                else egen long __lp_fa = group(__lp_stack __lp_T `strk')
                if "`mode'" != "one" egen long __lp_fb = group(__lp_stack __lp_Q `strk')
                if "`mode'" == "xt" xtset __lp_fa
            }
        }
        forvalues e = -`pre'/`post' {
            local d = `season' * (floor(`e' / `season') + 1)
            if `d' == 0 {
                matrix `R' = nullmat(`R') \ (0, `e', `e', 0, ., ., ., 0, 0, ., 1)
                matrix `BS' = nullmat(`BS') \ J(1, `ng', .)
                matrix `WS' = nullmat(`WS') \ J(1, `ng', .)
                continue
            }
            _lpddd_fit __lp_dy if __lp_rel == `e', `fitopt' `wopt'
            scalar `sb'  = r(b)
            scalar `sse' = r(se)
            scalar `sN'  = r(N)
            scalar `sdf' = r(df)
            scalar `sncl' = r(ncl)
            local fitrc = r(rc)
            matrix `tmpm' = r(bk)
            matrix `BS' = nullmat(`BS') \ `tmpm'
            matrix `tmpm' = r(wk)
            matrix `WS' = nullmat(`WS') \ `tmpm'
            if `fitrc' != 0 {
                local lastrc = `fitrc'
                local ++nfail
            }
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
                if `sncl' < . {
                    if `sncl' > `Nclmax' local Nclmax = `sncl'
                }
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
            if "`engine'" == "builtin" {
                if "`mode'" == "one" egen long __lp_fa = group(__lp_stack `strk')
                else if "`mode'" == "xt" egen long __lp_fa = group(__lp_stack __lp_g)
                else egen long __lp_fa = group(__lp_stack __lp_T `strk')
                if "`mode'" != "one" egen long __lp_fb = group(__lp_stack __lp_Q `strk')
                if "`mode'" == "xt" xtset __lp_fa
            }
        }
        forvalues j = 1/`nsp' {
            _lpddd_fit __lp_sp`j', `fitopt' `wopt'
            scalar `sb'  = r(b)
            scalar `sse' = r(se)
            scalar `sN'  = r(N)
            scalar `sdf' = r(df)
            scalar `sncl' = r(ncl)
            local fitrc = r(rc)
            matrix `tmpm' = r(bk)
            matrix `BS' = nullmat(`BS') \ `tmpm'
            matrix `tmpm' = r(wk)
            matrix `WS' = nullmat(`WS') \ `tmpm'
            if `fitrc' != 0 {
                local lastrc = `fitrc'
                local ++nfail
            }
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
                if `sncl' < . {
                    if `sncl' > `Nclmax' local Nclmax = `sncl'
                }
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
    display as text "Estimation: " as result cond("`engine'" == "reghdfe", "reghdfe", cond("`mode'" == "xt", "xtreg, fe", "regress"))
    display as text "Fixed effects: " as result `"`fe'"'
    display as text "Standard errors clustered by: " as result cond("`cluster'" == "", "`group'", "`cluster'")

    if "`table'" != "notable" {
        local nr = rowsof(`R')
        display _newline as text "{hline 14}{c TT}{hline 64}"
        display as text %13s "Periods" " {c |}" %10s "Estimate" %10s "Std. err." %8s "t" %7s "P>|t|" ///
            %21s "[`level'% conf. interval]" %8s "Obs."
        display as text "{hline 14}{c +}{hline 64}"
        forvalues i = 1/`nr' {
            if `R'[`i', 1] == 0 local lab = string(`R'[`i', 2])
            else local lab = string(`R'[`i', 2]) + " to " + string(`R'[`i', 3])
            if `R'[`i', 11] == 1 {
                display as text %13s "`lab'" " {c |}" as text %10s "0" as text "  (reference)"
            }
            else {
                display as text %13s "`lab'" " {c |}" as result %10.4f `R'[`i', 4] %10.4f `R'[`i', 5] ///
                    %8.2f `R'[`i', 6] %7.3f `R'[`i', 7] %10.4f `R'[`i', 8] %11.4f `R'[`i', 9] %8.0f `R'[`i', 10]
            }
            if `i' < `nr' {
                if `R'[`i', 1] == 0 & `R'[`i' + 1, 1] == 1 display as text "{hline 14}{c +}{hline 64}"
            }
        }
        display as text "{hline 14}{c BT}{hline 64}"
    }
    if `nfail' > 0 {
        display as error "Note: `nfail' regression(s) could not be estimated; last return code r(`lastrc')."
        if "`engine'" == "builtin" display as error "      With many fixed effects or weights that vary within a group, install reghdfe."
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
    if "`aggregate'" != "ols" {
        ereturn matrix stack_b = `BS'
        ereturn matrix stack_w = `WS'
    }
    if `Nclmax' > 0 ereturn scalar N_clust = `Nclmax'
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
    if "`weight'" != "" {
        ereturn local wtype "`weight'"
        ereturn local wexp "`exp'"
    }
    ereturn local engine "`engine'"
    ereturn local vce "cluster"
    ereturn local vcetype "Robust"
    if "`dd'" != "" ereturn local title "Stacked local-projection difference-in-differences"
    else ereturn local title "Stacked local-projection triple differences"
    ereturn local cmdline `"`cmdline'"'
    ereturn local cmd "lpddd"
end


* one regression: the coefficient on the treatment indicator, or an explicitly
* weighted average of stack-specific coefficients estimated in one regression
* (the average and its variance are computed from e(b) and e(V) in double precision)
program define _lpddd_fit, rclass
    version 15.1
    syntax varname [if] , d(varname) clv(varname) agg(string) nst(integer) ///
        engine(string) mode(string) [ fe(string) wtype(string) ]
    local wx ""
    if "`wtype'" != "" local wx "[`wtype' = __lp_w]"
    tempname b se N df ncl bk wk bf Vf tot est var
    scalar `b' = .
    scalar `se' = .
    scalar `N' = .
    scalar `df' = .
    scalar `ncl' = .
    matrix `bk' = J(1, `nst', .)
    matrix `wk' = J(1, `nst', .)

    if "`agg'" == "ols" local x "`d'"
    else {
        local x ""
        forvalues k = 1/`nst' {
            local x `x' __lp_D`k'
        }
    }
    if "`engine'" == "reghdfe" {
        capture reghdfe `varlist' `x' `if' `wx', absorb(`fe') vce(cluster `clv')
    }
    else if "`mode'" == "xt" {
        capture xtreg `varlist' `x' i.__lp_fb `if' `wx', fe vce(cluster `clv') nonest
    }
    else if "`mode'" == "two" {
        capture regress `varlist' `x' i.__lp_fa i.__lp_fb `if' `wx', vce(cluster `clv')
    }
    else {
        capture regress `varlist' `x' i.__lp_fa `if' `wx', vce(cluster `clv')
    }
    local rc = _rc

    if `rc' == 0 {
        scalar `N' = e(N)
        scalar `df' = e(df_r)
        scalar `ncl' = e(N_clust)
        if "`agg'" == "ols" {
            capture scalar `b' = _b[`d']
            if _rc == 0 scalar `se' = _se[`d']
            else scalar `b' = .
        }
        else {
            matrix `bf' = e(b)
            matrix `Vf' = e(V)
            local cn : colnames `bf'
            tempvar es
            quietly generate byte `es' = e(sample)
            scalar `tot' = 0
            forvalues k = 1/`nst' {
                matrix `wk'[1, `k'] = 0
                local nm : word `k' of `cn'
                local base = substr("`nm'", strrpos("`nm'", ".") + 1, .)
                if "`base'" == "__lp_D`k'" {
                    if `Vf'[`k', `k'] > 0 & `Vf'[`k', `k'] < . {
                        matrix `bk'[1, `k'] = `bf'[1, `k']
                        if "`agg'" == "equal" matrix `wk'[1, `k'] = 1
                        else {
                            quietly count if `es' & __lp_D`k' == 1
                            matrix `wk'[1, `k'] = r(N)
                        }
                        scalar `tot' = `tot' + `wk'[1, `k']
                    }
                }
            }
            if `tot' > 0 {
                scalar `est' = 0
                scalar `var' = 0
                forvalues k = 1/`nst' {
                    matrix `wk'[1, `k'] = `wk'[1, `k'] / `tot'
                }
                forvalues k = 1/`nst' {
                    if `wk'[1, `k'] > 0 {
                        scalar `est' = `est' + `wk'[1, `k'] * `bf'[1, `k']
                        forvalues l = 1/`nst' {
                            if `wk'[1, `l'] > 0 {
                                scalar `var' = `var' + `wk'[1, `k'] * `wk'[1, `l'] * `Vf'[`k', `l']
                            }
                        }
                    }
                }
                scalar `b' = `est'
                scalar `se' = sqrt(`var')
            }
        }
    }
    return scalar b = `b'
    return scalar se = `se'
    return scalar N = `N'
    return scalar df = `df'
    return scalar ncl = `ncl'
    return scalar rc = `rc'
    return matrix bk = `bk'
    return matrix wk = `wk'
end


program define _lpddd_parsesave, rclass
    version 15.1
    syntax anything(name=fn) [, replace]
    return local file `"`fn'"'
    return local replace "`replace'"
end
