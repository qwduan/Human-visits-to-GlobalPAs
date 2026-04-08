library(ggnewscale); library(DHARMa); library(readr) 
library(dplyr); library(mgcViz)
library(terra); library(sf)

fra_df <- read_csv('xxx/output/visit/home_fra_update.csv')
home <- fra_df %>% select(lat, lon, pop)

intensity_df <- read_csv('xxx/output/visit/df_inten_per_adj.csv') %>%
  left_join(home, by = c('visit_latitude'='lat', 'visit_longitude' = 'lon'))

sr_bird <- rast("xxx/Birds/Richness_10km.tif")
sr_mam <- rast("xxx/Mammals/Richness_10km.tif")
sr_amp <- rast("xxx/Amphibians/Richness_10km.tif")

sf_pts <- st_as_sf(intensity_df, coords = c("visit_longitude", "visit_latitude"), crs = 4326)

res <- 0.05
xmin <- min(intensity_df$visit_longitude, na.rm = TRUE)- 5*res/2
xmax <- max(intensity_df$visit_longitude, na.rm = TRUE)+ 5*res/2
ymin <- min(intensity_df$visit_latitude, na.rm = TRUE)- 5*res/2
ymax <- max(intensity_df$visit_latitude, na.rm = TRUE)+ 5*res/2

template_wgs <- rast(
  xmin = xmin, xmax = xmax, 
  ymin = ymin, ymax = ymax,
  resolution = c(res, res), 
  crs = "EPSG:4326" 
)

sf_pts$value_wgs <- intensity_df$sum_visit_time
r_wgs <- rasterize(vect(sf_pts), template_wgs, field = "value_wgs", fun = "mean", background = 0)

target_crs <- crs(sr_bird)
r_proj <- terra::project(r_wgs, target_crs, method = "bilinear")



adm0 <- st_read("xxx/World_Countries_Generalized.shp")
adm0 <- st_transform(adm0, crs = crs(r_wgs))
adm0_vect <- vect(adm0)

adm0_rast_orig <- rasterize(adm0_vect, r_wgs, field = "COUNTRY")
sum_orig <- zonal(r_wgs, adm0_rast_orig, fun = "sum", na.rm = TRUE)
colnames(sum_orig)[2] <- "visitsum_original"


adm0_proj <- st_transform(adm0, crs = crs(r_proj))
adm0_proj_vect <- vect(adm0_proj)
adm0_rast_proj <- rasterize(adm0_proj_vect, r_proj, field = "COUNTRY")
sum_proj <- zonal(r_proj, adm0_rast_proj, fun = "sum", na.rm = TRUE)
colnames(sum_proj)[2] <- "visitsum_projected"

country_fix <- adm0_proj %>%
  left_join(sum_orig, by = "COUNTRY") %>%
  left_join(sum_proj, by = "COUNTRY") %>%
  mutate(factor = visitsum_original / visitsum_projected)

country_fix_vect <- vect(country_fix)
zonal_factor <- rasterize(country_fix_vect, r_proj, field = "factor")

r_proj_corrected <- r_proj * zonal_factor

global(r_wgs, fun = 'sum', na.rm = TRUE)
global(r_proj, fun = 'sum', na.rm = TRUE)
global(r_proj_corrected , fun = "sum", na.rm = TRUE)


r_final <- resample(r_proj_corrected, sr_bird, method = "sum")
global(r_final , fun = "sum", na.rm = TRUE)


centroids <- as.points(r_final)
bird <- terra::extract(sr_bird, centroids)
mam <- terra::extract(sr_mam, centroids)
amp <- terra::extract(sr_amp, centroids)
visit <- terra::extract(r_final, centroids)


centroids_sf <- st_as_sf(centroids)
coords <- st_coordinates(centroids_sf)
centroids_df <- data.frame(
  lon = coords[, 1],
  lat = coords[, 2],
  sr_bird = bird[, 2],
  sr_mam = mam[, 2],
  sr_amp = amp[, 2],
  visit_intensity = visit[, 2]
)

nodata_map <- st_read("xxx/output/visit/map_country/country_non.shp")
nodata_map <- st_transform(nodata_map, crs = st_crs(centroids_sf))
idx_inside <- which(lengths(st_intersects(centroids_sf, nodata_map)) > 0)

centroids_df <- centroids_df[-idx_inside, ]

centroids_df$visit_intensity[is.na(centroids_df$visit_intensity)] <- 0


centroids_df <- centroids_df %>%
  mutate(sr = rowSums(across(c(sr_bird, sr_mam, sr_amp)), na.rm = TRUE))
centroids_df <- centroids_df %>%
  filter(!(is.na(sr_mam) & is.na(sr_bird) & is.na(sr_amp)))


com <- centroids_df
over <- read_csv('xxx/output/visit/factor10km_pa_intersect.csv')

over <- over %>%
  mutate(id = as.character(Id))

com$id <- as.character(com$id)
com_overlap <- com %>%
  left_join(
    over %>%
      group_by(id) %>%
      summarise(overlap = sum(fac10_over_area)),
    by = 'id'
  )

com_overlap$overlap[is.na(com_overlap$overlap)] <- 0
com_overlap$in_PA <- ifelse(com_overlap$overlap / 100 > 0.5, "Inside PA", "Outside PA")
com_overlap$in_PA_alt <- ifelse(com_overlap$overlap > 0, "Inside PA", "Outside PA")



hf_proj <- rast('xxx/output/visit/factor/hf.reproject.tif')
hf_final <- resample(hf_proj, sr_bird, method = "average")

slope_proj <- rast('xxx/output/visit/factor/slope.reproject.tif')
slope_final <- resample(slope_proj, sr_bird, method = "average")

ele_proj <- rast('xxx/output/visit/factor/ele.reproject.tif')
ele_final <- resample(ele_proj, sr_bird, method = "average")

tem_proj <- rast('xxx/output/visit/factor/tem.reproject.tif')
tem_final <- resample(tem_proj, sr_bird, method = "average")

ndvi_proj <- rast('xxx/output/visit/factor/ndvi.reproject.tif')
ndvi_final <- resample(ndvi_proj, sr_bird, method = "average")

land_proj <- rast('xxx/output/visit/factor/land.reproject.tif')
land_proj <- as.factor(land_proj)
land_final <- resample(land_proj, sr_bird, method = "mode")


pts <- vect(com_overlap, geom = c("lon", "lat"), crs = target_crs)
com_overlap$hf <- terra::extract(hf_final, pts)[,2]*0.01
com_overlap$slope <- terra::extract(slope_final, pts)[,2]
com_overlap$ele <- terra::extract(ele_final, pts)[,2]
com_overlap$tem <- terra::extract(tem_final, pts)[,2]
com_overlap$ndvi <- terra::extract(ndvi_final, pts)[,2]
com_overlap$land <- terra::extract(land_final, pts)[,2]

prec_orig <- rast("xxx/output/visit/factor/prec.tif") 
prec_reproj <- terra::project(prec_orig, target_crs, method = "bilinear")
zone_orig <- rasterize(adm0_vect, prec_orig, field = "COUNTRY")
zone_proj <- rasterize(adm0_proj_vect, prec_reproj, field = "COUNTRY")
sum_prec_orig <- zonal(prec_orig, zone_orig, fun = "sum", na.rm = TRUE)
sum_prec_proj <- zonal(prec_reproj, zone_proj, fun = "sum", na.rm = TRUE)
names(sum_prec_orig)[2] <- "prec_sum_orig"
names(sum_prec_proj)[2] <- "prec_sum_proj"
factor_df <- adm0_proj %>%
  left_join(sum_prec_orig, by = "COUNTRY") %>%
  left_join(sum_prec_proj, by = "COUNTRY") %>%
  mutate(prec_correction = prec_sum_orig / prec_sum_proj)

factor_raster <- rasterize(vect(factor_df), prec_reproj, field = "prec_correction")
prec_corrected <- prec_reproj * factor_raster

global(prec_orig, fun = 'sum', na.rm = TRUE)
global(prec_reproj, fun = 'sum', na.rm = TRUE)
global(prec_corrected, fun = "sum", na.rm = TRUE)

prec_fina <- resample(prec_corrected, sr_bird, method = "sum")
com_overlap$prec <- terra::extract(prec_fina, pts)[,2]




com <- com_overlap

library(mgcv)
vars <- c("sr", "logv", "ele", "slope", "tem", "prec", "ndvi", "hf", "lon", "lat")
com_clean <- com[complete.cases(com[, vars]), ]
data_in <- com_clean %>% filter(in_PA == 'Inside PA')

data_in$logsr <- log(data_in$sr)
data_in$land <- as.factor(data_in$land)
pos <- data_in$visit_intensity[!is.na(data_in$visit_intensity) & data_in$visit_intensity > 0]
c0 <- 0.5 * min(pos)
data_in$logv <- log(data_in$visit_intensity + c0)

fit <- gam(
  logsr ~ s(logv, k = 3, bs = "ts") +
    s(ele, k = 3, bs = "ts") +
    s(slope, k = 3, bs = "ts") +
    s(tem, k = 3, bs = "ts") +
    s(prec, k = 3, bs = "ts") +
    s(ndvi, k = 3, bs = "ts") +
    s(hf, k = 3, bs = "ts") +
    s(lon, lat, k = 3, bs = "ts"),
  data = data_in
)
summary(fit)
AIC(fit)
#gam.check(fit)  
concurvity(fit, full = TRUE)

sim_gam <- simulateResiduals(fittedModel = fit, plot = FALSE)
plot(sim_gam)
DHARMa::testResiduals(sim_gam)
DHARMa::testDispersion(sim_gam)


library(visreg)
v <- visreg(fit, "logv", scale = "response", partial = TRUE)
fit_df <- v$fit %>%
  mutate(
    # Use a small threshold to force the tiny numbers to 0
    visit_raw = exp(logv) - c0,
    
    sr_fit    = exp(visregFit),
    sr_lwr    = exp(visregLwr),
    sr_upr    = exp(visregUpr)
  )


ggplot(fit_df, aes(x = logv, y = sr_fit)) +
  geom_ribbon(aes(ymin = sr_lwr, ymax = sr_upr), fill = "#C8E6C9", alpha = 0.6) +
  geom_line(color = "#4CAF50", linewidth = 0.8) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )



##############

install.packages("gam.hp")  
library(gam.hp)
res_dev  <- gam.hp(fit, type = "dev")   
res_adj  <- gam.hp(fit, type = "adjR2") 

print(res_dev)
plot(res_dev)  


hp <- as.data.frame(res_dev$hierarchical.partitioning)
hp$term <- rownames(hp)

library(forcats)
# build plotting df
plot_df <- hp %>%
  transmute(
    term,
    individual = as.numeric(Individual),
    term_clean = gsub("^s\\(|,?\\s*k\\s*=.*\\)$", "", term),
    color = ifelse(grepl("\\blogv\\b", term), "#66B032", "#FF91A4")
  )
plot_df$term_clean <- factor(
  plot_df$term_clean,
  levels = plot_df$term_clean[order(plot_df$individual)]  # low→high (top will be high)
)
# bar plot (horizontal)
ggplot(plot_df, aes(x = individual, y = term_clean, fill = term_clean)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = setNames(plot_df$color, plot_df$term_clean)) +
  geom_text(aes(label = round(individual, 3)),
            hjust = -0.15, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Individual effect (explained deviance)",
       y = "Term") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

