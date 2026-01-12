setwd(getwd())

# Bring in needed libraries and functions
library("ggplot2")
source("Functions_Fo91.R")


# Bring in the csv file
olivine <- read.csv(
  "cat_olivine_demo.csv")
olivine <- subset(olivine, Fo < 91)
olivine <- subset(olivine, Fo > 82)
demo_alkalis <- 0.89
demo_silica <- 48.92
olivine$T_Fo91 <- NA

demo <- calculate_Fo91(olivine, 0.0001, 1, demo_alkalis, demo_silica)
demo_result <- demo[[1]]
demo_plot <- demo[[2]]
demo_plot


ggsave(filename = "Algorithm_stages.png",
       plot = demo_plot, width = 8, height = 3.5)

write.csv(demo_result, file = "calculated_Fo91_T.csv")

