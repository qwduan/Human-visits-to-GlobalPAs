library(sdmTMB)
library(INLA)

vars <- c("elevation", "slope", "hf", "road_den", "access", "build_den", "ndvi")

fishnet_pa <-  read_csv(file.path(folder, "fishnet_pa.csv"))


match <- read_csv('xxx/output/visit/covar_resample/match_var.csv')
pa_var <- fishnet_pa %>% left_join(match, by = 'fishnetid') %>%
  filter(!is.na(access))


#filter PAs where less than 95% of their area is covered by tiles that have data data, e.g. hf
pa_coverage_check <- pa_var %>%
  group_by(pa_id, pa_area) %>%
  summarise(
    total_overlap = sum(overlap_area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(overlap_ratio = total_overlap / pa_area) %>%
  filter(overlap_ratio >= 0.95) 
#0.985, 79073/80278

pa_var_filtered <- pa_var %>%
  semi_join(pa_coverage_check, by = "pa_id")


pa_summary <- pa_var_filtered %>%
  mutate(across(all_of(vars), ~ .x * overlap_area, .names = "w_{.col}")) %>%
  group_by(pa_id, pa_area, MIN_IUCN_rank, FIRST_ISO3) %>%
  summarise(
    total_overlap = sum(overlap_area, na.rm = TRUE),
    across(starts_with("w_"), ~ sum(.x, na.rm = TRUE) / total_overlap,
           .names = "{sub('w_', '', .col)}"),
    .groups = "drop"
  )

local_pop <- read_csv(file.path(folder, "local_100_pop.csv")) %>%
  dplyr::select(pa_id, SUM) %>%
  rename(local100pop = SUM)

pa_summary <- pa_summary %>% left_join(local_pop, by = 'pa_id')
pa_summary_f <- pa_summary %>%
  mutate(local100pop = replace_na(local100pop, 0))


gdp <- read_csv(file.path(folder, "gdp_2023.csv"))
hdi <- read_csv(file.path(folder, "hdi_2022.csv"))%>%
  mutate(iso3 = countrycode(Country, "country.name", "iso3c")) %>%
  select(-Country)
gover <- read_csv(file.path(folder, "government_2023.csv"))


pa_summary_f <- pa_summary_f %>% 
  left_join(gdp, by =c('FIRST_ISO3'= 'iso3')) %>%
  left_join(hdi, by =c('FIRST_ISO3'= 'iso3')) %>%
  left_join(gover, by =c('FIRST_ISO3'= 'iso3'))

pa <- read_csv(file.path(folder, "pa.csv")) %>%
  select(pa_id, centroid_lat,  centroid_lon)
pa_summary_f  <- pa_summary_f  %>% left_join(pa, by ='pa_id')


pa_summary_f$continent <- countrycode(pa_summary_f$FIRST_ISO3, origin = "iso3c", destination = "continent")
pa_summary_f$regions <- countrycode(pa_summary_f$FIRST_ISO3, origin = "iso3c", destination = "region")
pa_summary_f$region <- pa_summary_f$continent
pa_summary_f$region[
  pa_summary_f$continent == "Americas" & pa_summary_f$regions == "North America"
] <- "North America"

pa_summary_f$region[
  pa_summary_f$continent == "Americas" & pa_summary_f$regions == "Latin America & Caribbean"
] <- "Latin America & Caribbean"
table(pa_summary_f$region)

pa_summary_f <- pa_summary_f %>%
  mutate(
    HDI = as.numeric(HDI),
    ge = as.numeric(ge),
    cc = as.numeric(cc),
    rq = as.numeric(rq)
  )





visit_id <- read_csv("xxx/output/visit/v4/fishnet_humanvisit.csv") %>%
  dplyr::select(fishnetid, visit_latitude, visit_longitude)



######
visit_space <- read_csv('xxx/output/visit/v4/df_space.csv')
visit_index_space <- visit_space %>% left_join(visit_id, by = c('visit_latitude','visit_longitude')) %>% select(-...1)

for_space <- fishnet_pa %>%
  left_join(visit_index_space, by = "fishnetid")

visit_space_type <- for_space  %>%
  mutate(
    space_all = !is.na(category_code)
  ) %>%
  group_by(pa_id) %>%
  summarise(
    pa_area = first(pa_area),
    space_all = sum(overlap_area[space_all], na.rm = TRUE) / pa_area,
    .groups = "drop"
  )%>% 
  select(pa_id, space_all)


pa_summary_space <- pa_summary_f %>% 
  left_join(visit_space_type, by = 'pa_id')




pa_summary_space$iucn_group <- dplyr::case_when(
  pa_summary_space$MIN_IUCN_rank %in% c("1", "2", "3") ~ "Astrict",
  pa_summary_space$MIN_IUCN_rank %in% c("4", "5", "6", "7") ~ "Bless_strict",
  TRUE ~ "Cnot_classified"
)

predictors <- c("pa_area", "iucn_group", "elevation", "slope",
                "hf", "road_den", "access", "build_den", "ndvi", 
                "local100pop", "gdp",  "cc", "region")

data_space <- pa_summary_space %>%
  select(space_all, pa_id, centroid_lat,  centroid_lon, FIRST_ISO3, all_of(predictors)) %>%
  mutate(
    iucn_group = as.factor(iucn_group),
    FIRST_ISO3 = as.factor(FIRST_ISO3),
    region = as.factor(region)
  ) %>%
  na.omit()
write.csv(data_space, 'xxx/output/visit/v4/pa_modelvar_space.csv')




data_space <- read_csv('xxx/output/visit/v4/pa_modelvar_space.csv')
to_log <- c("pa_area", "access", "build_den", "road_den", "local100pop")
data_space <- data_space %>%
  mutate(across(
    all_of(c("pa_area", "access", "build_den", "road_den", "local100pop")),
    ~ log(.x + 0.1),
    .names = "log_{.col}"
  ))
space_data_scaled <- data_space %>%
  mutate(across(
    c(log_pa_area, elevation, slope, hf, log_road_den, log_access, log_build_den,
      ndvi, log_local100pop, gdp, cc),
    ~ as.numeric(scale(.x))
  ))
space_data_scaled$space_all[space_data_scaled$space_all > 1] <- 1


hist(space_data_scaled$space_all, breaks = 50, main = "Distribution of space_all", xlab = "space_all")




var_space <- c("log_pa_area", "iucn_group", "slope", "elevation","hf",
               "log_local100pop", "log_road_den", "log_access", "log_build_den", 
               "ndvi", "gdp",  "cc", "region")

form_space <- as.formula(paste("space_all ~", paste(var_space, collapse = " + "), "+ (1 | FIRST_ISO3)"))

space_data_scaled_sf <- st_as_sf(space_data_scaled, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
space_data_scaled_sf <- st_transform(
  space_data_scaled_sf,
  "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
)

coords2 <- st_coordinates(space_data_scaled_sf)
space_data_scaled$X_km <- coords2[, 1] / 1000
space_data_scaled$Y_km <- coords2[, 2] / 1000

mesh_space <- make_mesh(space_data_scaled, xy_cols = c("X_km", "Y_km"), cutoff = 50) 



space_data_scaled$y_zero <- ifelse(space_data_scaled$space_all == 0, 1, 0)
space_data_scaled$y_one <- ifelse(space_data_scaled$space_all == 1, 1, ifelse(space_data_scaled$space_all < 1 & space_data_scaled$space_all != 0, 0, NA))
space_data_scaled$y_proportion <- ifelse(space_data_scaled$space_all < 1 & space_data_scaled$space_all > 0, space_data_scaled$space_all, NA)


space_data_scaled$FIRST_ISO3 <- as.factor(space_data_scaled$FIRST_ISO3)
space_data_scaled$region <- as.factor(space_data_scaled$region)
space_data_scaled$iucn_group <- as.factor(space_data_scaled$iucn_group)


form_zero <- as.formula(paste("y_zero ~", paste(var_space, collapse = " + "), "+ (1 | FIRST_ISO3)"))

fit_zero <- sdmTMB(
  formula = form_zero,
  data = space_data_scaled,
  mesh = mesh_space,
  spatial = "on",
  family = binomial(link = "logit")
)


form_one <- as.formula(paste("y_one ~", paste(var_space, collapse = " + "), "+ (1 | FIRST_ISO3)"))
data_one <- subset(space_data_scaled, !is.na(y_one))

data_one_sf <- st_as_sf(data_one, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
data_one_sf <- st_transform(
  data_one_sf,
  "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
)
coords_one <- st_coordinates(data_one_sf)
data_one$X_km <- coords_one[, 1] / 1000
data_one$Y_km <- coords_one[, 2] / 1000
mesh_one <- make_mesh(data_one, xy_cols = c("X_km", "Y_km"), cutoff = 50)  

fit_one <- sdmTMB(
  formula = form_one,
  data = data_one,
  mesh = mesh_one,
  spatial = "on",
  family = binomial(link = "logit")
)



form_proportion <- as.formula(paste("y_proportion ~", paste(var_space, collapse = " + "), "+ (1 | FIRST_ISO3)"))
data_proportion <- subset(space_data_scaled, !is.na(y_proportion))

data_proportion_sf <- st_as_sf(data_proportion, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
data_proportion_sf <- st_transform(
  data_proportion_sf,
  "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
)
coords_proportion <- st_coordinates(data_proportion_sf)
data_proportion$X_km <- coords_proportion[, 1] / 1000
data_proportion$Y_km <- coords_proportion[, 2] / 1000
mesh_proportion <- make_mesh(data_proportion, xy_cols = c("X_km", "Y_km"), cutoff = 50)  


fit_proportion <- sdmTMB(
  formula = form_proportion,
  data = data_proportion,
  mesh = mesh_proportion,
  spatial = "on",
  family = Beta(link = "logit")
)


sim_zero <- simulate(fit_zero, nsim =100, type = "mle-mvn")
r_zero <- dharma_residuals(sim_zero, fit_zero, return_DHARMa = TRUE)
plot(r_zero)
#DHARMa::testResiduals(r_zero)
DHARMa::testDispersion(r_zero)
DHARMa::testZeroInflation(r_zero)

sim_one <- simulate(fit_one, nsim =100, type = "mle-mvn")
r_one <- dharma_residuals(sim_one, fit_one, return_DHARMa = TRUE)
plot(r_one)
#DHARMa::testResiduals(r_one)
DHARMa::testDispersion(r_one)
DHARMa::testZeroInflation(r_one)


sim_proportion <- simulate(fit_proportion, nsim =100, type = "mle-mvn")
r_proportion <- dharma_residuals(sim_proportion, fit_proportion, return_DHARMa = TRUE)
plot(r_proportion)
#DHARMa::testResiduals(r_sdm)
DHARMa::testDispersion(r_proportion)
DHARMa::testZeroInflation(r_proportion)





coef(fit_zero)
coef(fit_one)
coef(fit_proportion)
tidy(fit_proportion, "ran_pars")


coef_zero <- tidy(fit_zero, effects = "fixed",conf.int = TRUE)
coef_zero$model <- "Zero"
coef_one <- tidy(fit_one, effects = "fixed",conf.int = TRUE)
coef_one$model <- "One"
coef_prop <- tidy(fit_proportion, effects = "fixed",conf.int = TRUE)
coef_prop$model <- "Proportion"
coef_all <- rbind(coef_zero, coef_one, coef_prop)
write.csv(coef_all, 'xxx/output/visit/v4/space_sdm.csv')



############################

visit_intensity <- read_csv('xxx/output/visit/v4/df_inten_per.csv')
visit_index_inten <- visit_intensity %>% left_join(visit_id, by = c('visit_latitude','visit_longitude')) %>% select(-...1)

for_inten <- fishnet_pa %>%
  left_join(visit_index_inten, by = "fishnetid")

for_inten  <- for_inten  %>%
  mutate(
    sum_visit_time = replace_na(sum_visit_time, 0),
    rural_visit_time = replace_na(rural_visit_time, 0),
    sub_visit_time = replace_na(sub_visit_time, 0),
    urban_visit_time = replace_na(urban_visit_time, 0)
  )


for_inten <- for_inten  %>%
  mutate(
    inten_den = sum_visit_time / fishnet_area,
    w_inten_total = inten_den * overlap_area
  )

visit_inten_type <- for_inten %>%
  group_by(pa_id) %>%
  summarise(
    pa_area = first(pa_area),
    MIN_IUCN_rank = first(MIN_IUCN_rank),
    FIRST_ISO3 = first(FIRST_ISO3),
    total_overlap = sum(overlap_area, na.rm = TRUE),
    inten_all = sum(w_inten_total, na.rm = TRUE) / pa_area,
    .groups = "drop"
  ) %>% 
  select(pa_id, inten_all)




pa_summary_inten <- pa_summary_f %>% 
  left_join(visit_inten_type, by = 'pa_id')


pa_summary_inten$iucn_group <- dplyr::case_when(
  pa_summary_inten$MIN_IUCN_rank %in% c("1", "2", "3") ~ "Astrict",
  pa_summary_inten$MIN_IUCN_rank %in% c("4", "5", "6", "7") ~ "Bless_strict",
  TRUE ~ "Cnot_classified"
)


hist(pa_summary_inten$inten_all)


predictors <- c("pa_area", "iucn_group", "elevation", "slope",
                "hf", "road_den", "access", "build_den", "ndvi", 
                "local100pop", "gdp",  "cc", "region")

data_inten <- pa_summary_inten %>%
  select(inten_all, pa_id, centroid_lat,  centroid_lon, FIRST_ISO3, all_of(predictors)) %>%
  mutate(
    iucn_group = as.factor(iucn_group),
    FIRST_ISO3 = as.factor(FIRST_ISO3),
    region = as.factor(region)
  ) %>%
  na.omit()
write.csv(data_inten, 'xxx/output/visit/v4/pa_modelvar_inten.csv')
#77942



data_inten <- read_csv('xxx/output/visit/v4/pa_modelvar_inten.csv')
to_log <- c("pa_area", "access", "build_den", "road_den", "local100pop")
data_inten <- data_inten %>%
  mutate(across(
    all_of(c("pa_area", "access", "build_den", "road_den", "local100pop")),
    ~ log(.x + 0.1),
    .names = "log_{.col}"
  ))
inten_data_scaled <- data_inten %>%
  mutate(across(
    c(log_pa_area, elevation, slope, hf, log_road_den, log_access, log_build_den,
      ndvi, log_local100pop, gdp, cc),
    ~ as.numeric(scale(.x))
  ))



var_inten <- c("log_pa_area", "iucn_group", "slope", "elevation","hf",
          "log_local100pop", "log_road_den", "log_access", "log_build_den", 
          "ndvi", "gdp",  "cc", 'region')

form_inten <- as.formula(paste("inten_all ~", paste(var_inten, collapse = " + "), "+ (1 | FIRST_ISO3)"))


inten_scaled_sf <- st_as_sf(inten_data_scaled, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
inten_scaled_sf <- st_transform(
  inten_scaled_sf,
  "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
)

coords_inten <- st_coordinates(inten_scaled_sf)
inten_data_scaled$X_km <- coords_inten[, 1] / 1000
inten_data_scaled$Y_km <- coords_inten[, 2] / 1000

mesh_inten <- make_mesh(inten_data_scaled, xy_cols = c("X_km", "Y_km"), cutoff = 50)  


inten_data_scaled$FIRST_ISO3 <- as.factor(inten_data_scaled$FIRST_ISO3)
inten_data_scaled$region <- as.factor(inten_data_scaled$region)
inten_data_scaled$iucn_group <- as.factor(inten_data_scaled$iucn_group)

fit_inten <- sdmTMB(
  formula = form_inten,
  data = inten_data_scaled,
  mesh = mesh_inten,
  spatial = "on",
  family = tweedie()
)
sim <- simulate(fit_inten, nsim =100, type = "mle-mvn")
r_sdm <- dharma_residuals(sim, fit_inten, return_DHARMa = TRUE)
plot(r_sdm)
#DHARMa::testResiduals(r_sdm)
DHARMa::testDispersion(r_sdm)
DHARMa::testZeroInflation(r_sdm)


coef_inten <- tidy(fit_inten, conf.int = TRUE)

write.csv(coef_inten, 'xxx/output/visit/v4/inten_sdm.csv')
