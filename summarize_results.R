library(tidyverse)
setwd('~/Desktop/untracked/')
fs = dir(pattern='results')
d = tibble()
for(f in fs){
    d = read_csv(f, col_types = 'ccnc') %>%
        filter(! is.na(Year),
               ! is.na(Result)) %>%
        rename(region = Region, site = Site, year = Year, result = Result) %>%
        bind_rows(d)
}

successes = filter(d, result == 'Run Finished') %>% select(region, site, year)
failures = filter(d, result != 'Run Finished') %>% select(region, site, year)

write_csv(successes, 'new_models.csv')
