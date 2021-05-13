library(tidyverse)
library(StreamPULSE)
library(streamMetabolizer)

rating_cuve_ZQ = read.csv('~/git/streampulse/model/ZQ_data.csv')
pressure_sensor_offset_m = read.csv('~/git/streampulse/model/sensor_offsets.csv') %>%
    filter(site == 'NHC') %>%
    pull(offset_cm) %>%
    {. / 100}

nhc_data = StreamPULSE::request_data(
    sitecode = 'NC_NHC',
    startdate = '2016-01-01',
    enddate = '2021-05-12',
)

nhc_data_prepped = StreamPULSE::prep_metabolism(
    d = nhc_data,
    type = 'bayes',
    model = 'streamMetabolizer',
    interval = '15 min',
    zq_curve = list(sensor_height = pressure_sensor_offset_m,
                    Z = rating_cuve_ZQ[zq$site == 'NHC', 'level_m'],
                    Q = rating_cuve_ZQ[zq$site == 'NHC', 'discharge_cms'],
                    fit = 'power',
                    ignore_oob_Z = TRUE,
                    plot = TRUE),
    estimate_areal_depth = FALSE,
    estimate_PAR = TRUE,
    retrieve_air_pres = TRUE
)

nhc_q = as_tibble(nhc_data_prepped$data) %>%
    mutate(datetime_utc = streamMetabolizer::convert_solartime_to_UTC(
        any.solar.time = solar.time,
        longitude = -79.0018,
        time.type = 'mean solar'
    )) %>%
    select(datetime_utc, discharge_m3s = discharge)

write_csv(nhc_q, '~/Desktop/nhc_Q.csv')
