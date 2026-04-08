#####################################################estimate visit space and visit intensity within PAs
#fishnet_pa: overlap between pixels and PAs

library(dplyr); library(readr); library(ggplot2); library(tidyr); library(countrycode)

folder <- "xxx/output/visit/covar_resample"
pa <- read_csv(file.path(folder, "pa_greater1.csv"))

#filter urban pixels
fishnet_pa <- read_csv(file.path(folder, "fishnet_pa_up.csv")) %>%
  mutate(fishnetid = fishernet_fishnetid, 
         fishnet_area = fishernet_fishnet_area)

fishnet_pa$pa_ratio <- fishnet_pa$overlap_area/fishnet_pa$fishernet_fishnet_area

fishnet_urban <- read_csv(file.path(folder, "fishent_urban.csv")) %>%
  mutate(
    fishnetid   = fishernet_fishnetid,
    fishnet_area = fishernet_fishnet_area
  ) %>%
  group_by(fishnetid) %>%
  summarise(
    sum_urban    = sum(overlap_area, na.rm = TRUE),
    fishnet_area = first(fishnet_area),
    urban_area   = sum_urban / fishnet_area,
    .groups = "drop"
  ) %>%
  select(fishnetid, urban_area)

fishnet_pa  <- fishnet_pa %>% filter(overlap_area > 0) %>%
  left_join(fishnet_urban, by = 'fishnetid') %>%
  mutate(pa_ratio = overlap_area / fishnet_area)

fishnet_pa <- fishnet_pa %>%
  filter(is.na(urban_area) | urban_area < 0.25)

visit_id <- read_csv("xxx/output/visit/fishnet_humanvisit.csv") %>%
  mutate(fishnetid = fishernet_fishnetid) %>%
  dplyr::select(fishnetid, visit_latitude, visit_longitude)

df <- fishnet_pa  %>%
  left_join(visit_id, by = "fishnetid")

df <- df %>%
  mutate(size = cut(
    pa_area,
    breaks = c(1, 10, 20, 30, 40, 50, 100, 1000, 100000, Inf),
    labels = c(
      "1–10", "10–20", "20–30", "30–40", "40–50",
      "50–100", "100–1,000", "1,000–100,000", ">100,000"
    ),
    right = FALSE
  ))

df$continent <- countrycode(df$FIRST_ISO3, origin = "iso3c", destination = "continent")
table(df$continent)
df$continent[df$FIRST_ISO3 == 'CCK'] <- "Oceania" 

df$regions <- countrycode(df$FIRST_ISO3, origin = "iso3c", destination = "region")

df$region <- df$continent
df$region[
  df$continent == "Americas" & df$regions == "North America"
] <- "North America"

df$region[
  df$continent == "Americas" & df$regions == "Latin America & Caribbean"
] <- "Latin America & Caribbean"
table(df$region)

#install.packages("WDI")
library(WDI)
# Get country income levels from the World Bank
income<- WDI(indicator = "NY.GDP.PCAP.CD", extra = TRUE)

country_income <- income %>% 
  filter(year == '2024')%>%
  dplyr::select(income, iso3c)
country_income$income[country_income$iso3c == "CZE"] <- "High income"
country_income$income[country_income$iso3c == "GUF"] <- "High income"
country_income$income[country_income$iso3c == 'VNM'] <- 'Lower middle income'

df <- df %>%
  left_join(country_income, by = c('FIRST_ISO3'='iso3c'))
df$income <- factor(df$income, levels = c("High income", "Upper middle income", 
                                                "Lower middle income", "Low income"))


visit_inten <- read_csv('xxx/output/visit/df_inten_per_adj.csv')
visit_space <- read_csv('xxx1/output/visit/df_space.csv')

visit_index_inten <- df %>% left_join(visit_inten, by = c('visit_latitude','visit_longitude')) %>% select(-...1)
visit_index_space <- df %>% left_join(visit_space, by = c('visit_latitude','visit_longitude')) %>% select(-...1)

overlap_path <- "xxx/output/figures/2.overlap_pa"

visit_index_inten <- visit_index_inten %>%
  mutate(
    sum_visit_time = replace_na(sum_visit_time, 0),
    urban_visit_time = replace_na(urban_visit_time, 0),
    sub_visit_time = replace_na(sub_visit_time, 0),
    rural_visit_time = replace_na(rural_visit_time, 0)
  ) 

visit_index_space <- visit_index_space %>%
  mutate(
    rural_flag = if_else(category_code %in% c("1", "4", "5", "7"), 1L, 0L),
    sub_flag   = if_else(category_code %in% c("2", "4", "6", "7"), 1L, 0L),
    urban_flag = if_else(category_code %in% c("3", "5", "6", "7"), 1L, 0L),
    all_flag   = if_else(!is.na(category_code), 1L, 0L),
    
    all_overlap   = all_flag   * overlap_area,
    rural_overlap = rural_flag * overlap_area,
    sub_overlap   = sub_flag   * overlap_area,
    urban_overlap = urban_flag * overlap_area
  ) 

pa_info <- visit_index_space %>%
  group_by(pa_id) %>%
  summarise(
    area_pa = first(pa_area),
    region = first(region),
    income = first(income),
    size = first(size),
    IUCN_rank = first(MIN_IUCN_rank),
    .groups = "drop"
  )

calc_space_ratio <- function(data, pa_info, col_overlap) {
   data %>%
    group_by(pa_id) %>%
    summarise(
      space = sum(.data[[col_overlap]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
  right_join(pa_info, by = "pa_id") %>%
    mutate(
      space = replace_na(space, 0)
    )
}

all_ext <- calc_space_ratio(visit_index_space, pa_info, "all_overlap")
rural_ext <- calc_space_ratio(visit_index_space, pa_info, "rural_overlap")
sub_ext <- calc_space_ratio(visit_index_space, pa_info, "sub_overlap")
urban_ext <- calc_space_ratio(visit_index_space, pa_info, "urban_overlap")

all_ext %>%
  summarise(
    total_overlap = sum(space, na.rm = TRUE),
    total_pa = sum(area_pa, na.rm = TRUE),
    ratio = total_overlap / total_pa,
    .groups = "drop"
  )

summarise_overlap <- function(data, group_var) {
  data %>%
    group_by({{ group_var}}) %>%
    summarise(
      total_overlap = sum(space, na.rm = TRUE),
      total_pa = sum(area_pa, na.rm = TRUE),
      ratio = total_overlap / total_pa,
      .groups = "drop"
    )
}

rank_ex_rural <- summarise_overlap(rural_ext, IUCN_rank)%>%
  mutate(urban = "rural")
rank_ex_sub <- summarise_overlap(sub_ext, IUCN_rank)%>%
  mutate(urban = "sub")
rank_ex_urban <- summarise_overlap(urban_ext, IUCN_rank)%>%
  mutate(urban = "urban")
rank_ex_all <- summarise_overlap(all_ext, IUCN_rank)

summary_ranking <- bind_rows(rank_ex_rural, rank_ex_sub, rank_ex_urban)

ggplot(summary_ranking , aes(x = as.factor(IUCN_rank), y = ratio, fill = urban)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(
    rural = '#31625c',
    sub = '#fcd168',
    urban = '#b92336'
  ))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(overlap_path, "space_rank.pdf"),
  plot = last_plot(), 
  width = 8,
  height = 7
)

region_ext_rural <- summarise_overlap(rural_ext, region)%>%
  mutate(urban = "rural")
region_ext_sub <- summarise_overlap(sub_ext, region)%>%
  mutate(urban = "sub")
region_ext_urban <- summarise_overlap(urban_ext, region)%>%
  mutate(urban = "urban")
region_ex_all <- summarise_overlap(all_ext,  region)

space_region <- bind_rows(region_ext_rural, region_ext_sub, region_ext_urban)

global_space <- space_region %>%
  group_by(urban) %>%
  summarise(total_overlap = sum(total_overlap, na.rm = TRUE), 
            total_pa = sum(total_pa, na.rm=TRUE),
            ratio = total_overlap/total_pa,.groups = "drop") %>%
  mutate(region = "Global")


# Combine global summary with original region summaries
space_region_all <- bind_rows(space_region, global_space)
space_region_all <- space_region_all %>%
  mutate(region = factor(region, levels = c("Global", setdiff(unique(region), "Global"))))

ggplot(space_region_all, aes(x = region, y = ratio, fill = urban)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(
    rural = '#31625c',
    sub = '#fcd168',
    urban = '#b92336'
  ))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = file.path(overlap_path, "space_region.pdf"),
  plot = last_plot(), 
  width = 8,
  height = 7
)

size_ext_rural <- summarise_overlap(rural_ext, size)%>%
  mutate(urban = "rural")
size_ext_sub <- summarise_overlap(sub_ext, size)%>%
  mutate(urban = "sub")
size_ext_urban <- summarise_overlap(urban_ext, size)%>%
  mutate(urban = "urban")

summary_size <- bind_rows(size_ext_rural, size_ext_sub, size_ext_urban)

ggplot(summary_size, aes(x = size, y = ratio, fill = urban)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(
    rural = '#31625c',
    sub = '#fcd168',
    urban = '#b92336'
  ))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = file.path(overlap_path, "space_size.pdf"),
  plot = last_plot(), 
  width = 8,
  height = 7
)

income_ext_rural <- summarise_overlap(rural_ext, income)%>%
  mutate(urban = "rural")
income_ext_sub <- summarise_overlap(sub_ext, income)%>%
  mutate(urban = "sub")
income_ext_urban <- summarise_overlap(urban_ext, income)%>%
  mutate(urban = "urban")
income_ex_all <- summarise_overlap(all_ext,  income)
summary_income <- bind_rows(income_ext_rural, income_ext_sub, income_ext_urban)%>% filter(!is.na(income))

ggplot(summary_income, aes(x = as.factor(income), y = ratio, fill = urban)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(
    rural = '#31625c',
    sub = '#fcd168',
    urban = '#b92336'
  ))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = file.path(overlap_path, "space_income.pdf"),
  plot = last_plot(), 
  width = 8,
  height = 7
)


######################################intensity
visit_index_inten <- visit_index_inten %>%
  mutate(
    inten_den_rural = rural_visit_time / fishnet_area,
    inten_den_sub = sub_visit_time/ fishnet_area,
    inten_den_urban = urban_visit_time / fishnet_area,
    inten_den_all = sum_visit_time / fishnet_area,
    
    inten_over_rural = inten_den_rural * overlap_area,
    inten_over_sub = inten_den_sub * overlap_area,
    inten_over_urban = inten_den_urban * overlap_area,
    inten_over_all = inten_den_all * overlap_area
  )

df_inten_all <- visit_index_inten %>%
  group_by(pa_id) %>%
  summarise(
    inten_overlap = sum(inten_over_all, na.rm = TRUE),
    .groups = "drop"
  )

df_inten_all <- pa_info %>%
  left_join(df_inten_all, by = "pa_id") %>%
  mutate(
    inten_overlap = replace_na(inten_overlap, 0),
    intensity_ratio = inten_overlap / area_pa
  )

df_inten_all %>% 
  summarise(
    total_intensity = sum(inten_overlap, na.rm = TRUE),
    total_area = sum(area_pa, na.rm = TRUE),
    ratio = total_intensity / total_area,
    .groups = "drop"
  )

df_inten_long <- visit_index_inten %>%
  select(pa_id, pa_area, region, income, size, MIN_IUCN_rank,
         inten_over_rural, inten_over_sub, inten_over_urban) %>%
  pivot_longer(
    cols = starts_with("inten_over_"),
    names_to = "urban",
    values_to = "inten_over"
  ) %>%
  mutate(
    urban = case_when(
      urban == "inten_over_rural" ~ "rural",
      urban == "inten_over_sub" ~ "sub",
      urban == "inten_over_urban" ~ "urban"
    )
  )

calc_intensity_overlap <- function(data, pa_info) {
  intensity <- data %>%
    group_by(pa_id, urban) %>%
    summarise(
      inten_overlap = sum(inten_over, na.rm = TRUE),
      .groups = "drop"
    )
  
  pa_info %>%
    crossing(urban = c("rural", "sub", "urban")) %>%
    left_join(intensity, by = c("pa_id", "urban")) %>%
    mutate(
      inten_overlap = replace_na(inten_overlap, 0),
      intensity_ratio = inten_overlap / area_pa
    )
}

intensity_df <- calc_intensity_overlap(df_inten_long, pa_info)

summarise_intensity <- function(data, group_var) {
  data %>%
    group_by({{ group_var}}, urban) %>%
    summarise(
      total_intensity = sum(inten_overlap, na.rm = TRUE),
      total_pa = sum(area_pa, na.rm = TRUE),
      ratio = total_intensity / total_pa,
      .groups = "drop"
    )
}

rank <- summarise_intensity(intensity_df,  IUCN_rank)

ggplot(rank, aes(x = as.factor(IUCN_rank), y = ratio, fill = urban)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(
    rural = '#31625c',
    sub = '#fcd168',
    urban = '#b92336'
  ))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = file.path(overlap_path, "inten_rank.pdf"),
  plot = last_plot(), 
  width = 8,
  height = 7
)
rank_summary <- df_inten_all %>%
  group_by(IUCN_rank) %>%
  summarise(
    total_intensity = sum(inten_overlap, na.rm = TRUE),
    total_area = sum(area_pa, na.rm = TRUE),
    ratio = total_intensity / total_area,
    .groups = "drop"
  )


region <- summarise_intensity(intensity_df, region)

inten_global_summary <- region %>%
  group_by(urban) %>%
  summarise(total_intensity = sum(total_intensity, na.rm = TRUE), 
            total_pa = sum(total_pa, na.rm=TRUE),
            ratio = total_intensity/total_pa,.groups = "drop") %>%
  mutate(region = "Global")


inten_region_all <- bind_rows(region, inten_global_summary)
inten_region_all <- inten_region_all %>%
  mutate(region = factor(region, levels = c("Global", setdiff(unique(region), "Global"))))
ggplot(inten_region_all, aes(x = region, y = ratio, fill = urban)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(
    rural = '#31625c',
    sub = '#fcd168',
    urban = '#b92336'
  ))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = file.path(overlap_path, "inten_region.pdf"),
  plot = last_plot(), 
  width = 8,
  height = 7
)

region_summary <- df_inten_all %>%
  group_by(region) %>%
  summarise(
    total_intensity = sum(inten_overlap, na.rm = TRUE),
    total_area = sum(area_pa, na.rm = TRUE),
    ratio = total_intensity / total_area,
    .groups = "drop"
  )


size <- summarise_intensity(intensity_df,  size)

ggplot(size, aes(x = size, y = ratio, fill = urban)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(
    rural = '#31625c',
    sub = '#fcd168',
    urban = '#b92336'
  ))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = file.path(overlap_path, "inten_size.pdf"),
  plot = last_plot(), 
  width = 8,
  height = 7
)

income <-summarise_intensity(intensity_df, income)
income <- income %>% filter(!is.na(income))
ggplot(income, aes(x = as.factor(income), y = ratio, fill = urban)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(
    rural = '#31625c',
    sub = '#fcd168',
    urban = '#b92336'
  ))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = file.path(overlap_path, "inten_income.pdf"),
  plot = last_plot(), 
  width = 8,
  height = 7
)
income_summary <- df_inten_all %>%
  group_by(income) %>%
  summarise(
    total_intensity = sum(inten_overlap, na.rm = TRUE),
    total_area = sum(area_pa, na.rm = TRUE),
    ratio = total_intensity / total_area,
    .groups = "drop"
  )

