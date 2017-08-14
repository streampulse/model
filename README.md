# Model
Streampulse metabolism modeling code.

change.

# Workflow
Fitting metabolism models with streampulse follows two steps:

1. Downloading/formatting data from the streampulse platform
2. Fitting the chosen model
  * `model_name = "streamMetabolizer"` runs [streamMetabolizer](https://github.com/USGS-R/streamMetabolizer)
    - `model_type="bayes"` fits a Bayesian model (recommended)
    - `model_type="mle"` fits a maximum likelihood model
  * `model_name = "BASE"` runs [BASE](https://github.com/dgiling/BASE)

The code for these functions is contained in three files: `sp_functions.R`, `gapfill_functions.R`, and `BASE_functions.R` (which are copied from the `dgiling` repo).

# Code
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

# Model type for streamMetabolizer
# We recommend the Bayesian model, but you can also fit "mle", which runs much faster.
model_type <- "bayes"
# Which modeling framework to use
# "streamMetabolizer" is default, can also use "BASE"
model_name <- "streamMetabolizer"

# Download data from streampulse and prepare for metabolism modeling
fitdata <- prep_metabolism(sitecode = "NC_Eno",
    startdate = "2016-01-01", enddate = "2017-01-01",
    type = model_type, model = model_name, fillgaps = TRUE)

# Fit metabolism model
modelfit <- fit_metabolism(fitdata)

# Gather metabolism predictions
predictions <- predict_metabolism(modelfit)
```

# Contributing
We welcome outside contributions as pull requests that team members can review.

Project members can follow these steps:
1. Create a personal branch on the streampulse GitHub page.
2. Clone the repo. In the terminal, if you have ssh set up: `$ git clone git@github.com:streampulse/model.git`
3. Checkout your branch: `$ git checkout <branch-name>`
    and rebase to bring your branch up to date with master before you start editing the files: `$ git rebase master`
4. Make your edits, then push back to your branch:
		`$ git commit -am “<my message>”`
    `$ git push origin <branch-name>`
5. On GitHub, create a pull request to bring your changes back into the master branch (and maybe ask someone else to look at them).
