# Projective Averages for Summarizing Redistricting Ensembles

A recurring challenge in the application of redistricting simulation algorithms
lies in extracting useful summaries and comparisons from a large ensemble of
districting plans. Researchers often compute summary statistics for each
district in a plan, and then study their distribution across the plans in the
ensemble. This approach discards rich geographic information that is inherent
in districting plans. We introduce the projective average, an operation that
projects a district-level summary statistic back to the underlying geography
and then averages this statistic across plans in the ensemble. Compared to
traditional district-level summaries, projective averages are a powerful tool
for geographically granular, sub-district analysis of districting plans along a
variety of dimensions.  However, care must be taken to account for variation 
within redistricting ensembles, to avoid misleading conclusions.  We propose
and validate a multiple-testing procedure to control the probability of
incorrectly identifying outlier plans or regions when using projective
averages.

## Replication

To replicate the figures and analyses in the paper, run the scripts in `replication/` in order:

``` r
lapply(sort(Sys.glob("replication/*.R")), source)
```

Then run `quarto render paper/proj-avg.qmd` to generate the paper.

## Software

The methods described in the paper are implemented in the [**redist**](https://alarm-redist.org/redist/) software.
As of the time of writing (i.e., before version 4.3 is released), you will need the development version installed.
