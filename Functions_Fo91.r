# Partition coefficient Mg between olivine and liquid
D_ol_liq_mg <- function(P, H2O, T_c, alkali) {
  
  phase1 <- -2.158 + 55.09 * (P/T_c) - (0.06213) * H2O + (4430/T_c) + 
    (0.05115) * alkali
  
  phase2 <- exp(phase1)
  
  return(phase2)
  
}

# Partition coefficient Fe between olivine and liquid
D_ol_liq_fe <- function(P, H2O, T_c, alkali, SiO2) {
  
  phase1 <- -3.300 + 47.57 * (P/T_c) - (0.05192) * H2O + (3344/T_c) + 
    (0.05595) * alkali + (0.01633) * SiO2
  
  phase2 <- exp(phase1)
  
  return(phase2)
  
}

# Cation fraction of Mg in the liquid
X_liq_mg <- function(x , y) {
  
  result <- x / y
  
  return(result)
  
}

# Cation fraction of Fe in the liquid
X_liq_fe <- function(x , y) {
  
  result <- x / y
  
  return(result)
  
}

# Olivine composition in equilibrium with the liquid
X_ol_mg_div_X_ol_fe <- function(x, y) {
  
  mg <- (x/(x + y))*100
  result <-(100*mg)/(0.30*(100-mg)+mg)
  return(result)
}

ol_mg1 <- function(x, y, z) {
  
  result <- (z*x)/(0.3*y)
  return(result)
}

# Solve for partition coefficient of iron
D_for_fe <- function(x, y, z) {
  
  result <- (0.667 - x*y) / z
  return(result)
}

#Solve for T with partition coefficient of magnesium
new_T_mg <- function(D, P, H2O, alkali) {
  
  result <- (55.09*P+4430)/
    (log(D)+2.158+0.06213*H2O-0.05115*alkali)
  return(result)
}

#Solve for T with partition coefficient of iron
new_T_fe <- function(D, P, H2O, alkali, SiO2) {
  
  result <- (47.57*P+3344)/
    (log(D)+3.3+0.05192*H2O-0.05595*alkali-0.01633*SiO2)
  return(result)
}


calculate_Fo91 <- function(Loc, pressure, water, alkalis, silica) {
  
  temps <- c(0)
  fors <- c(0)
  samples <- c(0)
  pla <- c(0)
  plot_df_counter = 1
  df_for_plotting <- data.frame(Fo=double(),
                                T_C=double(),
                                Grain=character())
  
  for (i in 1:length(Loc$Study)) {
    
    id <- Loc$Grain[i]
    forsterite <- Loc$Fo[i]
    temperature <- Loc$T_C[i]
    pl <- Loc$Placement[i]
    ol_mg2 <- Loc$MgO[i]
    ol_fe2 <- Loc$FeO[i]
    sample <- Loc[i,][c(12:14, 16:23, 25)]
    olivine_other_elements <- sum(sample, na.rm = T)
    ol_mg <- ol_mg2/(olivine_other_elements+ol_fe2+ol_mg2)
    ol_fe <- ol_fe2/(olivine_other_elements+ol_fe2+ol_mg2)
    olivine_other_elements <- olivine_other_elements/(olivine_other_elements+ol_fe2+ol_mg2)
    df_for_plotting[plot_df_counter, 1] <- forsterite
    df_for_plotting[plot_df_counter, 2] <- temperature
    df_for_plotting[plot_df_counter, 3] <- id
    plot_df_counter <- plot_df_counter + 1
    
    
    while (forsterite < 91) {
      # Partition coefficient from the observed temperature
      D_mg <- D_ol_liq_mg(pressure, water, temperature, alkalis)
      D_fe <- D_ol_liq_fe(pressure, water, temperature, alkalis, silica)
      
      # Based on the partition coefficient the cation proportions
      # of the liquid phase are calculated
      liq_mg <- X_liq_mg(ol_mg, D_mg)
      liq_fe <- X_liq_fe(ol_fe, D_fe)
      other_elements <- 1-liq_mg-liq_fe
      
      # Forsterite in equilibrium with this melt can be calculated with this
      forsterite <- X_ol_mg_div_X_ol_fe(liq_mg, liq_fe)
      df_for_plotting[plot_df_counter, 1] <- forsterite
      df_for_plotting[plot_df_counter, 3] <- id
      
      # The cation proportions of the olivine in equilibrium with the new melt
      # are calculated here
      ol_mg <- ol_mg1(liq_mg, liq_fe, ol_fe)
      ol_fe <- liq_fe*D_fe
      
      # Simple mixing calculation can be done as such
      liq_mg2 <- (0.99*liq_mg+0.01*ol_mg)
      liq_fe2 <- (0.99*liq_fe+0.01*ol_fe)
      liq_mg <- liq_mg2/(other_elements+liq_fe2+liq_mg2)
      liq_fe <- liq_fe2/(other_elements+liq_fe2+liq_mg2)
      other_elements <- other_elements/(other_elements+liq_fe2+liq_mg2)
      
      # New partition coefficients are calculated with the cation proportions
      D_mg <- ol_mg/liq_mg
      D_fe <- (ol_fe/liq_fe)
      
      temperature <- new_T_mg(D_mg, pressure, water, alkalis)
      T2 <- new_T_fe(D_fe, pressure, water, alkalis, silica)
      
      # Simultaneously solving the partition coefficient for the ol/liq
      # and the olivine stoichiometry.
      while(round(D_mg*liq_mg+D_fe*liq_fe,2) != round(0.667, 2)){
        
        if(D_mg*liq_mg+D_fe*liq_fe < 0.667) {
          D_mg <- D_mg + 0.001
        } else {
          
          D_mg <- D_mg - 0.0001
          temperature <- new_T_mg(D_mg, pressure, water, alkalis)
          T2 <- new_T_fe(D_fe, pressure, water, alkalis, silica)
          #print(round(D_mg*liq_mg+D_fe*liq_fe,2))
          #print(c(temperature, T2))
          while (round(temperature, 2) != round(T2, 2)) {
            if (temperature > T2) {
              D_fe <- D_fe - 0.00001
              T2 <- new_T_fe(D_fe, pressure, water, alkalis, silica)
            } else if(T2 > temperature) {
              D_mg <- D_mg - 0.00001
              temperature <- new_T_mg(D_mg, pressure, water, alkalis)
            }
          }
        }
      }
      
      ol_mg <- liq_mg*D_mg
      ol_fe <- liq_fe*D_fe
      
      df_for_plotting[plot_df_counter, 2] <- temperature
      plot_df_counter <- plot_df_counter + 1
      
    }
    
    temps <- append(temps, round(temperature))
    fors <- append(fors, round(forsterite,0))
    samples <- append(samples, id)
    pla <- append(pla, pl)
    Loc$T_Fo91[i] <- round(temperature)
    
  }
  temps <- temps[-1]
  fors <- fors[-1]
  samples <- samples[-1]
  pla <- pla[-1]
  result <- data.frame("Temperatures" = temps,
                       "Fo" = fors,
                       "Grain" = samples,
                       "Placement" = pla)
  #print(df_for_plotting)
  
  # Creating a plot to visualise the algorithm stages
  plot_stages <- ggplot(data = df_for_plotting, aes(x = Fo,
                                                    y = T_C,
                                                    colour = Grain)) +
    geom_point(size = 2) +
    geom_path(linewidth = 1) +
    scale_color_discrete() +
    ylab(bquote(T[cryst]*"")) +
    xlab("Fo") +
    theme_classic()
  plot_stages
  
  return(list(Loc, plot_stages))
}




