# Streampulse metabolism modeling code

## Metabolism modeling workflow
Fitting metabolism models with streampulse follows three steps:

1. Downloading/formatting data from the streampulse platform (with `prep_metabolism()`)
2. Fitting the chosen model (with `fit_metabolism()`)
  * `model_name = "streamMetabolizer"` runs [streamMetabolizer](https://github.com/USGS-R/streamMetabolizer)
    - `model_type="bayes"` fits a Bayesian model (recommended)
    - `model_type="mle"` fits a maximum likelihood model
  * `model_name = "BASE"` runs [BASE](https://github.com/dgiling/BASE)
3. Predicting metabolic rates (with `predict_metabolism()`)

The code for these functions is contained in three files: `sp_functions.R`, `gapfill_functions.R`, and `BASE_functions.R` (which is adapted from the `dgiling/BASE` repo).

## Only getting data from streampulse?
If you just want to get data from streampulse into R, you can use `sp_data(sitecode, startdate, enddate, variables, flags, token)`
* `sitecode` is a site name, like "regionID_siteID"
* `startdate` and `enddate` are YYYY-MM-DD strings, e.g., "1983-12-09"
* `variables` is a vector of c("variable_one", ..., "variable_n")
* `flags` is logical, to include flag data (default = FALSE)
* `token` is a private data token (default = NULL), only relevant if you are accessing embargoed data (none of the core streampulse data are embargoed)

This function queries the RESTful API

## Code
This example is also available as code in [`model_streampulse.R`](https://github.com/streampulse/model/blob/master/model_streampulse.R).

```r
# You need to have an internet connection for this code to function.
# REQUIRED packages
library(coda)
library(dplyr)
library(httr)
library(jsonlite)
library(R2jags)
library(streamMetabolizer)
library(tidyr)

# REQUIRED code
# These source our StreamPULSE functions from GitHub
# In the future, turn these codes into a package and just source the package...
sp_functions <- GET("https://raw.githubusercontent.com/streampulse/model/master/sp_functions.R")
eval(parse(text = content(sp_functions, as="text", encoding="UTF-8")), envir= .GlobalEnv)

# Which modeling framework to use
# "streamMetabolizer" is default, can also use "BASE"
model_name <- "streamMetabolizer"

# Model type for streamMetabolizer
# We recommend the Bayesian model, but you can also fit "mle", which runs much faster.
model_type <- "bayes"

# Download data from streampulse and prepare for metabolism modeling
fitdata <- prep_metabolism(sitecode = "NC_Eno",
    startdate = "2016-01-01", enddate = "2017-01-01",
    type = model_type, model = model_name, fillgaps = TRUE)

# Fit metabolism model
modelfit <- fit_metabolism(fitdata)

# Gather metabolism predictions
predictions <- predict_metabolism(modelfit)
```

## Contributing
We welcome outside contributions as pull requests that team members can review.

Project members can follow these steps:
1. Create a personal branch on the streampulse GitHub page.
2. Clone the repo. In the terminal, if you have ssh set up: `$ git clone git@github.com:streampulse/model.git`
3. Before you start editing:
    * Checkout your branch `$ git checkout <branch-name>`
    * Rebase to bring your branch up to date with master: `$ git rebase origin master`
4. Make your edits, then push back to your branch:
		* `$ git commit -am “<my message>”`
    * `$ git push origin <branch-name>`
5. On GitHub, create a pull request to bring your changes back into the master branch (and maybe ask someone else to look at them).
