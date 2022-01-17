# GAM/GLMM doodles and scraps

## Creating a complete model using covariates such as bin lag and timeseries richness ####

# Create a giant data frame for individual communities

model_input <- frame_builder_function(Fish_Communities_B, 1, "B")

# This is a GLMM with site or timeseriesID treated as a random intercept
# 

covariate_full_glmm <- glmer(novel ~  log(bin.lag) + 
                               log(TimeSeries_Length) + 
                               (1|site) + 
                               log(position) + 
                               log(richness),
                             data=model_input, family=binomial)

# If the random intercept does not explain any variance, run as a GLM with fixed effects only instead.

covariate_full_glmm <- glm(novel ~  bin.lag + 
                             TimeSeries_Length + 
                             position + 
                             richness,
                           data=model_input, family=binomial)


# Remove non-significant factors
covariate_full_glmm <- glm(novel ~  log(bin.lag) + 
                             log(position) + 
                             log(richness),
                           data=model_input, family=binomial)


# Run a GAM with ybp as smoothing terms

model_input$ybp <- as.numeric(2021 - model_input$bins)

covariate_gam <- gam(novel ~ log(bin.lag) + 
                       log(position) + 
                       log(richness) +
                       s(ybp, bs="tp", k= -1),
                     data=model_input, family=betar(), method = "REML")


estimated_values <- predict_gam(covariate_gam)

estimated_values %>%
  ggplot(aes(ybp, plogis(fit))) +
  geom_smooth_ci(ybp)


summary(covariate_full_glmm)
summary(covariate_gam)

check <- gam(novel ~ bins, data = model_input)


a <- by(data = model_input$novel, INDICES = model_input$ybp, FUN = mean)
ybp <- as.vector(names(a))
new <- as.data.frame(cbind(as.vector(ybp), as.numeric(a)))
colnames(new) <- c("ybp", "prob")

plot(new$ybp, new$prob)
lines(loess.smooth(new$prob, new$ybp))

lines(t ,g$fitted.values,col=3)



saveRDS(model_input_B_1, "./outputs/GLM_input_B_1.rds")


