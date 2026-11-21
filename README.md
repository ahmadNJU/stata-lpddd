# lpddd

Stacked local-projection triple differences for Stata, with a seasonal base period and estimates for spans of periods.

Version 0.1.0 (22 November 2026). This is a first version for testing.

## What it does

`lpddd` estimates the dynamic effect of a treatment that starts at different dates for different groups and affects only some units within a treated group. Examples: a tariff on some products (groups) from one country (eligible units) but not from other countries; a state policy that applies to one category of persons.

It combines:

- **Stacking.** Each treated cohort is put in a stack with clean control groups; time is counted from the cohort's first treated period. A group treated earlier is never a control for one treated later.
- **Local projections.** One regression per period relative to treatment, on the change of the outcome from a base period.
- **Triple difference.** The regressor is treated group x eligible unit; fixed effects remove what is common to a group and to eligible or ineligible units, by stack.
- **Seasonal base period.** With `season(12)` every month is compared with the same calendar month of the last year before treatment.
- **Spans.** With `spans()` the outcome is added up over a span of periods (for example months 12 to 23) and compared with the same seasons before treatment.

A difference-in-differences is available with the option `dd`.

## Install

Requires Stata 14 or later and `reghdfe`:

```stata
ssc install reghdfe
ssc install ftools
```

From GitHub:

```stata
net install lpddd, from("https://raw.githubusercontent.com/ahmadNJU/stata-lpddd/main") replace
net get lpddd, from("https://raw.githubusercontent.com/ahmadNJU/stata-lpddd/main")
help lpddd
```

`net get` copies `lpddd_example.do` to the current folder.

## Example

```stata
* monthly imports by product and supplier; the tariff applies to supplier China
lpddd imports, unit(prod_supp) time(month) group(product) cohort(first_tariff) ///
    eligible(china) log season(12) pre(36) post(35) graph

* spans of months, comparisons within products and sectors
lpddd imports, unit(prod_supp) time(month) group(product) cohort(first_tariff) ///
    eligible(china) log season(12) noevent spans(-12 -1 0 2 3 11 12 23 24 35) ///
    gfe strata(chapter)
```

## Test

```stata
do lpddd_example.do
```

The do-file simulates data with a known effect and checks `lpddd` against two estimates made by hand: the triple difference of four cell means for one stack and one month, and a stacked regression with group and sector effects for a span. It stops with an error if they differ. If `lpdid` is installed it also displays `lpddd, dd` beside `lpdid, nevertreated`.

## Files

| File | Purpose |
|---|---|
| `lpddd.ado` | the command |
| `lpddd.sthlp` | help file |
| `lpddd.pkg`, `stata.toc` | package description for `net install` |
| `lpddd_example.do` | simulated example and tests |

## Planned

- response per unit of a treatment intensity (for example a tariff rate)
- contributions of unit types to the group total and the replacement rate
- covariates and reweighting of control groups

## References

- Cengiz, D., Dube, A., Lindner, A., and Zipperer, B. (2019). The effect of minimum wages on low-wage jobs. *Quarterly Journal of Economics*, 134, 1405-1454.
- Dube, A., Girardi, D., Jordà, Ò., and Taylor, A. M. (2025). A local projections approach to difference-in-differences. *Journal of Applied Econometrics*, 40(7). https://doi.org/10.1002/jae.70000
- Hsieh, M. H. (2026). Stacked triple differences. arXiv:2604.22982.
- Jordà, Ò. (2005). Estimation and inference of impulse responses by local projections. *American Economic Review*, 95, 161-182.
- Wing, C., Freedman, S. M., and Hollingsworth, A. (2024). Stacked difference-in-differences. NBER Working Paper 32054.

## Author

Ahmad Nawaz, Business School, Nanjing University; University of Sahiwal, Pakistan ahmadnawaz@uosahiwal.edu.pk
