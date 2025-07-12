# READ ALL DATA
# Hengxing Zou


########## Initialization ##########


library(tidyverse)
library(patchwork)

library(foreach)
library(doSNOW)

cl = makeCluster(4, outfile = "")
registerDoSNOW(cl)

color_2 = c("#2972B3", "#D98E04")
color_3 = c("#464DAE", "#CC7252", "#850394")
color_4 = c("#2972B3", "#734002", "#D98E04", "#99C8F2")

# Change the directory according to your working environment

database_dir = "Database/"
generated_dir = "GeneratedData/"


########## Read All Data ##########


# Read community-weighted metrics of all species

All_Weighted_Metrics = read_csv(paste0(generated_dir, "CWMetrics_Filtered.csv")) %>% 
  rename("Hand.Wing.Index" = "Hand-Wing.Index")

# Read environmental data

Bioclim = read_csv(paste0(database_dir, "Bioclim.csv"))

# Read anthrome data

Anthromes = read_csv(paste0(database_dir, "Anthromes_Grid_Years_Prop.csv")) %>% 
  pivot_wider(names_from = Anthrome, values_from = Coverage) %>% 
  rename(year = Year) %>% 
  # bin data into larger land use categories
  mutate(Settlements = `11` + `12` + `22` + `23` + `24`, 
         Agriculture = `31` + `32` + `33` + `34` + `41` + `42` + `43`, 
         Cultured = `51` + `52` + `53` + `54`, 
         Wild = `61` + `62` + `63`) %>% 
  select(-(7:25), -LONGD, -LATD, -AREA_SQ_KM, -Grid_ID)


########## Process CWM Data ##########


# Join CW Metrics with climate and anthrome data

cwm = All_Weighted_Metrics %>% 
  inner_join(Bioclim, by = c("Region" = "ST_12", "year" = "year")) %>% 
  inner_join(Anthromes, by = c("Region" = "ST_12", "year" = "year"))

cwm_mean = cwm %>% 
  filter(Metric == "CWM")

# Select traits

md_traits = colnames(All_Weighted_Metrics)[c(12, 14, 17, 16, 18:19, 21)] # for adding corrected generation lengths
md_traits_log = c("PC1_beak", "PC1_wing", "relative_bill_length", "relative_wing_length",
                  "log_Mass", "log_clutch_size", "corr_GenLength", "relative_gen_length")

# Trait labels

trait_labels = c(`PC1_beak` = "Beak PC1", 
                 `PC1_wing` = "Wing PC1", 
                 `relative_bill_length` = "Relative Beak Length",
                 `relative_wing_length` = "Relative Wing Length", 
                 `Mass` = "Mass", 
                 `litter_or_clutch_size_n` = "Clutch Size", 
                 `GenLength` = "Generation Length", 
                 `log_Mass` = "log(Mass)", # after transformation
                 `log_clutch_size` = "log(Clutch Size)", # after transformation
                 `log_GenLength` = "log(Generation Length)", # after transformation
                 `corr_GenLength` = "Corrected Generation Length",
                 `relative_gen_length` = "Relative Generation Length"
                )
