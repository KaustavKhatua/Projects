
# Loading data
data(Grunfeld, package = "AER")

# Converting the data in panel data
library(plm)
panel_data <- pdata.frame(Grunfeld, index = c("firm", "year"))


# Pooled Model
pooled_model <- lm(invest ~ value + capital, data = Grunfeld)
summary(pooled_model)


# Plotting dependent variable against dependent variables.
library(ggplot2)
library(grid)
library(gridExtra)

value_graph <- ggplot(Grunfeld[c(1:60, 201:220), ], aes(x = value, y = invest, color = factor(firm))) +
  geom_line() + theme(legend.position = "none", plot.margin = unit(c(1, 0, 1, 1), "lines")) + 
  labs(y = "Invest", x = "Value")

capital_graph <- ggplot(Grunfeld[c(1:60, 201:220), ], aes(x = capital, y = invest, color = factor(firm))) + 
  geom_line() + theme(legend.direction = "horizontal", plot.margin = unit(c(1, 0, 1, 1), "lines")) + 
  labs(y = "Invest", x = "Capital")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(capital_graph)

capital_graph <- capital_graph + theme(legend.position = "none")

grid.arrange(value_graph, capital_graph, legend, ncol = 2, nrow = 2, layout_matrix =
               rbind(c(1, 2), c(3, 3)), widths = c(2.7, 2.7), heights = c(2.5, 0.2), top = 
               textGrob("Invest as Function of Regressors", gp = gpar(fontsize = 14, font = 1)))


# Durbin Watson Test for Pooled Model.
library(plm)
pdwtest(invest ~ value + capital, data = panel_data, method = "pooling")


# Fixed Effects Model
fixed_effects_model <- lm(invest ~ value + capital + firm - 1, data = Grunfeld)  # - 1 to fit model without intercept.
summary(fixed_effects_model)

# Durbin-Watson test for Fixed Effects Model
pdwtest(invest ~ value + capital, data = panel_data, model = "within")


# Random Effects Model
random_effects_model <- plm(invest ~ value + capital, data = panel_data, model = "random")
summary(random_effects_model)

# R-Squared of Random Effects Model using different formula.
r <- random_effects_model$residuals
1 - (sum(r^2) / sum(  (Grunfeld[, 1] - ( sum(Grunfeld[, 1]) / nrow(Grunfeld)) )^2  ))


# Durbin-Watson Test for Random Effects Model
pdwtest(invest ~ value + capital, data = panel_data, model = "random")



# Pooled vs Random Effects Model
((sum(pooled_model$residuals^2) - sum(fixed_effects_model$residuals^2)) / 
    sum(fixed_effects_model$residuals^2)) * ((220 - 2 - 11) / 10)



# Fixed vs Random Effects Model (Hausman Test)
# First we have to fit Fixed Effects Model using plm package, as phtest function takes models fitted
# using plm function.
p_fixed_effects_model <- plm(invest ~ value + capital, data = panel_data, model = "within")
phtest(p_fixed_effects_model, random_effects_model)



# Handling Autocorrelation
# Prais-Winsten transformation.
rho_vector <- 1:2

for(i in 1:11){
  subset_data <- Grunfeld[((i - 1) * 20 + 1):((i - 1) * 20 + 20), ]
  e <- lm(invest ~ value + capital - 1, data = subset_data)$residuals
  numerator <- 0
  for(j in 2:length(e)){
    numerator <- numerator + (e[j] * e[j - 1])
  }
  rho_vector[i] <- numerator / sum(e^2)
}

transformed_data <- Grunfeld

for(i in 1:11){
  rho <- rho_vector[i]
  for(j in ((i - 1) * 20 + 20):((i - 1) * 20 + 2)){
    transformed_data[j, 1:3] <- transformed_data[j, 1:3] - rho * transformed_data[j - 1, 1:3]
    # transformed_data[j, 1:3] <- transformed_data[j, 1:3] * (1 - rho^2)^0.5
  }
  transformed_data[(i - 1) * 20 + 1, 1:3] <- transformed_data[(i - 1) * 20 + 1, 1:3] * (1 - rho^2)^0.5
}

transformed_panel_data <- pdata.frame(transformed_data, index = c("firm", "year"))

# Fitting model on transformed data.
transformed_random_effects_model <- plm(invest ~ value + capital, data = transformed_panel_data, 
                                        model = "random")
summary(transformed_random_effects_model)


# Durbin-Watson Test on the transformed model.
pdwtest(invest ~ value + capital, data = transformed_panel_data, model = "random")