{smcl}
{* *! version 0.1.1  22nov2026}{...}
{viewerjumpto "Syntax" "lpddd##syntax"}{...}
{viewerjumpto "Description" "lpddd##description"}{...}
{viewerjumpto "Options" "lpddd##options"}{...}
{viewerjumpto "Method" "lpddd##method"}{...}
{viewerjumpto "Examples" "lpddd##examples"}{...}
{viewerjumpto "Stored results" "lpddd##results"}{...}
{viewerjumpto "References" "lpddd##references"}{...}
{viewerjumpto "Author" "lpddd##author"}{...}
{title:Title}

{phang}
{bf:lpddd} {hline 2} Stacked local-projection triple differences with a seasonal base period


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:lpddd} {depvar} {ifin} {weight}{cmd:,}
{opth u:nit(varname)} {opth t:ime(varname)} {opth gro:up(varname)} {opth coh:ort(varname)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth u:nit(varname)}}panel unit, observed over time{p_end}
{synopt:{opth t:ime(varname)}}time period (integers, for example a monthly date){p_end}
{synopt:{opth gro:up(varname)}}group at which treatment is assigned; one or more units per group{p_end}
{synopt:{opth coh:ort(varname)}}first treated period of the group; 0 or missing if never treated{p_end}

{syntab:Design}
{synopt:{opth el:igible(varname)}}0/1 indicator of the units that are treated when their group is treated; gives a triple difference{p_end}
{synopt:{opt dd}}difference-in-differences instead (treated against control groups); not with {opt eligible()}{p_end}
{synopt:{opt cont:rol(never|notyet)}}control groups: never treated (default), or never treated and not yet treated within the window{p_end}

{syntab:Time}
{synopt:{opt pre(#)}}number of periods before treatment in the period-by-period estimates; default 0{p_end}
{synopt:{opt post(#)}}number of periods after the first treated period; default 0{p_end}
{synopt:{opt seas:on(#)}}length of the seasonal cycle (12 for monthly data, 4 for quarterly); default 1, no season{p_end}
{synopt:{opth span:s(numlist)}}pairs of first and last period of spans to be estimated as one outcome{p_end}
{synopt:{opt noev:ent}}do not estimate period by period (only {opt spans()}){p_end}

{syntab:Outcome}
{synopt:{opt log}}use the natural log of {depvar}; for spans, the log of the total over the span{p_end}
{synopt:{opt ihs}}use the inverse hyperbolic sine instead, which keeps zeros{p_end}
{synopt:{opt fill:zero}}in spans, count periods without an observation as zero{p_end}

{syntab:Model}
{synopt:{opt gfe}}stack-by-group fixed effects in place of stack-by-treated-status effects{p_end}
{synopt:{opth str:ata(varname)}}interact the fixed effects with this variable (for example a sector){p_end}
{synopt:{opt abs:orb(absvars)}}further fixed effects, passed to {cmd:reghdfe}{p_end}
{synopt:{opt agg:regate(ols|cohort|equal)}}how the stacks are combined; default {cmd:ols}{p_end}
{synopt:{opth cl:uster(varname)}}cluster variable; default is {opt group()}{p_end}
{synopt:{opt eng:ine(reghdfe|builtin)}}estimation command; default is {cmd:reghdfe} if installed, otherwise built-in commands{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt gra:ph}}plot the period-by-period estimates with their confidence interval{p_end}
{synopt:{opt sav:ing(filename[, replace])}}save the results as a Stata dataset{p_end}
{synopt:{opt notab:le}}suppress the table{p_end}
{synoptline}
{p 4 6 2}
{cmd:aweight}s, {cmd:pweight}s and {cmd:fweight}s are allowed and must be constant within a unit.{p_end}
{p 4 6 2}
{cmd:lpddd} runs without other packages. If {cmd:reghdfe} (Correia) is installed it is used, which is faster with many fixed effects; {opt absorb()} requires it.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:lpddd} estimates the dynamic effect of a treatment that starts at different dates for
different groups and that affects only some units within a treated group. Examples are a
tariff imposed on some products (groups) from one country (eligible units) but not from
other countries, or a state policy (group) that applies to one category of persons only.

{pstd}
The command combines four elements.

{phang2}1. {it:Stacking.} For each treated cohort (all groups first treated in the same
period) a stack is formed that holds this cohort and clean control groups. Time is counted
from the first treated period of the cohort. A group treated earlier is never a control
for a group treated later.

{phang2}2. {it:Local projections.} One regression is run for each period relative to
treatment. The outcome is the change between that period and a base period before
treatment.

{phang2}3. {it:Triple difference.} The regressor is treated group x eligible unit. The
fixed effects remove what is common to a group and what is common to eligible or
ineligible units, separately in each stack.

{phang2}4. {it:Seasonal base period.} With {opt season(s)} each period is compared with
the same season of the last cycle before treatment, so seasonal movements cancel in the
difference. The last cycle before treatment (periods -s to -1) is the reference and is
zero by construction. With {cmd:season(1)} the base is the last period before treatment,
as in the usual event study.

{pstd}
With {opt spans()} the command also estimates effects for spans of periods (for example
months 12 to 23). The outcome is first added up over the span and then compared with the
same seasons before treatment. Spans are less noisy than single periods, and with
{opt log} a unit stays in the sample as long as its total over the span is positive.

{pstd}
With {opt dd} the third difference is dropped and the command estimates a stacked
local-projection difference-in-differences for the units in the sample.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth unit(varname)} is the panel unit. In a triple difference a group holds at least one
eligible and one ineligible unit, for example a product imported from the taxed country
and the same product imported from all other countries.

{phang}
{opth time(varname)} is the period. It must contain integers with one unit per period, as
Stata monthly, quarterly or yearly dates do. {opt unit()} and {opt time()} must identify
the observations.

{phang}
{opth group(varname)} is the level at which the treatment date is set. Standard errors are
clustered by this variable unless {opt cluster()} is given.

{phang}
{opth cohort(varname)} is the first treated period of the group, in the units of
{opt time()}. It is 0 or missing for groups that are never treated, and it must be
constant within a group.

{dlgtab:Design}

{phang}
{opth eligible(varname)} marks with 1 the units that receive the treatment when their
group is treated, and with 0 the other units. It must be constant within a unit.

{phang}
{opt dd} estimates a difference-in-differences between treated and control groups. Use
{it:if} to choose the units, for example the eligible units alone or the ineligible units
alone. Estimating the ineligible units alone shows whether they responded to the
treatment of the eligible units.

{phang}
{opt control(never)} uses the never-treated groups as controls in every stack.
{opt control(notyet)} adds groups that are first treated after the end of the stack's
window.

{dlgtab:Time}

{phang}
{opt pre(#)} and {opt post(#)} set the periods of the period-by-period estimates, from
-{it:pre} to +{it:post}, where 0 is the first treated period. Periods before -{it:s} are
compared with the same season of the reference cycle; a coefficient different from zero
there indicates different movements before treatment.

{phang}
{opt season(#)} is the length of the seasonal cycle. A period {it:e} is compared with
period {it:e} - {it:s}(floor({it:e}/{it:s}) + 1), which lies between -{it:s} and -1.

{phang}
{opth spans(numlist)} gives pairs of numbers, the first and last period of each span, for
example {cmd:spans(-12 -1 0 2 3 11 12 23)}. A span after treatment is compared with the
same seasons of the reference cycle: with {cmd:season(12)}, months 3 to 11 are compared
with months -9 to -1 and months 12 to 23 with months -12 to -1. A span before treatment is
compared with the same periods one cycle earlier (or as many full cycles as the span
needs): months -12 to -1 are compared with months -24 to -13. A span may not contain
periods both before and after treatment. Note the difference in sign for periods before
treatment: period by period, an earlier period is measured against the reference cycle;
in spans, a span is measured against the periods before it.

{phang}
{opt noevent} skips the period-by-period estimates.

{dlgtab:Outcome}

{phang}
{opt log} and {opt ihs} transform the outcome. Period by period, the outcome is
f(y) in the period less f(y) in the base period. In spans it is f(average of y over the
span) less f(average of y over the reference periods). Without either option the outcome
is used as it is, and a variable that is already in logs can be used directly for the
period-by-period estimates.

{phang}
{opt fillzero} applies to spans. By default a unit enters a span only if it is observed
in every period of the span and of its reference periods. With {opt fillzero} missing
periods count as zero, which suits data in which a missing record means no transaction.

{dlgtab:Model}

{pstd}
By default the regression for a triple difference contains stack-by-treated-status and
stack-by-eligibility fixed effects. Together with the regressor these saturate the four
cells of each stack (treated or control group, eligible or ineligible unit).

{phang}
{opt gfe} replaces the stack-by-treated-status effects by stack-by-group effects, so that
eligible and ineligible units are compared within the same group.

{phang}
{opth strata(varname)} interacts the fixed effects with a variable such as a sector, so
that comparisons are made within its categories. With {opt gfe} only the
stack-by-eligibility effects are interacted.

{phang}
{opt absorb(absvars)} adds fixed effects in the notation of {cmd:reghdfe} and requires
{cmd:reghdfe}. For spans the variables should not vary over time within a unit.

{phang}
{opt aggregate(ols)} estimates one coefficient from the pooled stacks. It is a weighted
average of the stack-specific effects with positive weights that depend on the size and
composition of the stacks. {opt aggregate(cohort)} estimates one coefficient per stack in
a single regression and averages them with weights equal to the number of treated units
of each stack in the estimation sample. {opt aggregate(equal)} gives every stack the same
weight. The standard error of the average comes from the same regression and accounts
for control groups that appear in several stacks.

{phang}
{opth cluster(varname)} sets the cluster variable. Control groups appear in several
stacks; clustering by {opt group()} or a coarser variable accounts for this.

{phang}
{opt engine(reghdfe)} estimates with {cmd:reghdfe}; {opt engine(builtin)} uses
{cmd:regress}, or {cmd:xtreg, fe} with {opt gfe}. The coefficients are the same. Standard
errors can differ slightly because the commands count degrees of freedom differently.
With {opt engine(builtin)} the fixed effects other than the group effects enter as
indicator variables, so many strata need a large {cmd:matsize} in Stata 15, and with
{opt gfe} weights must be constant within a group. A note is shown if a regression
could not be estimated.

{dlgtab:Reporting}

{phang}
{opt graph} draws the period-by-period estimates.
{opt saving(filename[, replace])} saves the results (variables type, from, to, b, se, t,
p, ll, ul, N, ref).


{marker method}{...}
{title:Method}

{pstd}
For stack {it:g}, unit {it:i} and period {it:e} relative to treatment, the regression is

{p 8 8 2}
f(y[i,e]) - f(y[i,base(e)]) = tau[e] x T[i] x Q[i] + fixed effects + error

{pstd}
where T marks the treated groups of the stack and Q the eligible units. One regression is
run for each {it:e}, on all stacks together. Period by period this is the local-projection
form of the event study (Jord{c a'} 2005; Dube, Girardi, Jord{c a'} and Taylor 2025) applied
to stacked data (Cengiz, Dube, Lindner and Zipperer 2019). With the default fixed effects
it is the saturated stacked triple-difference regression of Hsieh (2026), whose pooled
coefficient is a positively weighted average of the stack-specific triple differences;
{opt aggregate(cohort)} and {opt aggregate(equal)} are the explicit weights discussed
there. The seasonal base period and the spans are additions of this command.

{pstd}
The estimates rest on the assumption that, without treatment, the difference between
eligible and ineligible units would have moved in the same way in treated and control
groups, and that ineligible units and control groups are not affected by the treatment.


{marker examples}{...}
{title:Examples}

{pstd}Monthly imports by product and supplier; the tariff applies to supplier China{p_end}
{phang2}{cmd:. lpddd imports, unit(prod_supp) time(month) group(product) cohort(first_tariff) eligible(china) log season(12) pre(36) post(35) graph}{p_end}

{pstd}Spans of months, products weighted, comparisons within sectors and within products{p_end}
{phang2}{cmd:. lpddd imports [aweight = imports2017], unit(prod_supp) time(month) group(product) cohort(first_tariff) eligible(china) log season(12) noevent spans(-24 -13 -12 -1 0 2 3 11 12 23 24 35) gfe strata(chapter)}{p_end}

{pstd}Other suppliers alone: did they gain in the taxed products?{p_end}
{phang2}{cmd:. lpddd imports if china == 0, unit(prod_supp) time(month) group(product) cohort(first_tariff) dd log season(12) noevent spans(12 23 24 35) strata(chapter)}{p_end}

{pstd}Every tariff list with the same weight{p_end}
{phang2}{cmd:. lpddd imports, unit(prod_supp) time(month) group(product) cohort(first_tariff) eligible(china) log season(12) post(35) aggregate(equal)}{p_end}

{pstd}A simulated example is in {cmd:lpddd_example.do} (type {cmd:net get lpddd}).{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:lpddd} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}largest number of observations in a regression{p_end}
{synopt:{cmd:e(N_clust)}}largest number of clusters in a regression{p_end}
{synopt:{cmd:e(n_stacks)}}number of stacks (treated cohorts){p_end}
{synopt:{cmd:e(season)}}length of the seasonal cycle{p_end}
{synopt:{cmd:e(pre)}}periods before treatment{p_end}
{synopt:{cmd:e(post)}}periods after treatment{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:lpddd}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of the outcome{p_end}
{synopt:{cmd:e(design)}}{cmd:ddd} or {cmd:dd}{p_end}
{synopt:{cmd:e(cohorts)}}first treated periods of the stacks{p_end}
{synopt:{cmd:e(control)}}{cmd:never} or {cmd:notyet}{p_end}
{synopt:{cmd:e(aggregate)}}{cmd:ols}, {cmd:cohort} or {cmd:equal}{p_end}
{synopt:{cmd:e(transform)}}{cmd:none}, {cmd:log} or {cmd:ihs}{p_end}
{synopt:{cmd:e(absorb)}}fixed effects of the regressions{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable{p_end}
{synopt:{cmd:e(vce)}}{cmd:cluster}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust}{p_end}
{synopt:{cmd:e(engine)}}{cmd:reghdfe} or {cmd:builtin}{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(title)}}title of the estimation{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}estimates; names m{it:#} (period -#), p{it:#} (period +#), span{it:#}{p_end}
{synopt:{cmd:e(V)}}variances on the diagonal; the estimates come from separate regressions, so covariances are not available{p_end}
{synopt:{cmd:e(results)}}one row per period or span: type, from, to, b, se, t, p, ll, ul, N, ref{p_end}
{synopt:{cmd:e(stack_b)}}with {opt aggregate(cohort)} or {opt aggregate(equal)}: estimate of each stack, one row per row of {cmd:e(results)}{p_end}
{synopt:{cmd:e(stack_w)}}the weights given to the stacks{p_end}

{pstd}
{cmd:e(sample)} is not set: the estimates come from many regressions on the stacked data, not on the data in memory.


{marker references}{...}
{title:References}

{phang}
Cengiz, D., A. Dube, A. Lindner, and B. Zipperer. 2019. The effect of minimum wages on
low-wage jobs. {it:Quarterly Journal of Economics} 134: 1405-1454.

{phang}
Correia, S. 2017. Linear models with high-dimensional fixed effects: An efficient and
feasible estimator. Working paper.

{phang}
Dube, A., D. Girardi, {c O'}. Jord{c a'}, and A. M. Taylor. 2025. A local projections
approach to difference-in-differences. {it:Journal of Applied Econometrics} 40(7).

{phang}
Hsieh, M. H. 2026. Stacked triple differences. arXiv:2604.22982.

{phang}
Jord{c a'}, {c O'}. 2005. Estimation and inference of impulse responses by local
projections. {it:American Economic Review} 95: 161-182.

{phang}
Wing, C., S. M. Freedman, and A. Hollingsworth. 2024. Stacked difference-in-differences.
NBER Working Paper 32054.


{marker author}{...}
{title:Author}

{pstd} Ahmad Nawaz {p_end}
{pstd} Business School, Nanjing University, China{p_end}
{pstd} Department of Economics, University of Sahiwal, Pakistan {p_end}
{pstd} ahmadnawaz@uosahiwal.edu.pk {p_end}


{title:Also see}

{psee}
{helpb reghdfe}, {helpb lpdid}, {helpb stackdid} (if installed)
{p_end}
