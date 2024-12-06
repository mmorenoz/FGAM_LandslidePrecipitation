# The code illustrates the procedures described in the manuscript--
# Functional regression for space-time prediction of precipitation-induced shallow landslides in South Tyrol, Italy
# by authors Mateo Moreno, Luigi Lombardo, Stefan Steger, Lotte de Vugt,
# Thomas Zieher, Alice Crespi, Francesco Marra, Cees van Westen and Thomas Opitz


# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install package renv for reproductibility
# install.packages("renv")

# install packages with respective version
renv::restore()

# call packages
list.packages <- c("mgcv", "refund" ,"sf", "tidyverse", "mgcViz", "forcats")
vapply(list.packages, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)
remove(list.packages)

# ggplot theme
theme_plot <- function() {
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        plot.title = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,3,1,1), "lines"),
        axis.ticks.length.x = unit(.3, "cm"))
}


# LOADING DATA ------------------------------------------------------------

# load cum precipitation
d <- readRDS("./DATA/d_5_0cumforward_5day_t0123.Rds") |> 
  tidyr::drop_na(day1.00:day5.23)
table(d$bin)

# precipitation preparation
d_prec <- d |> 
  sf::st_drop_geometry() |> 
  dplyr::select(day1.00:day5.23) |> 
  `colnames<-`(round(seq(1,6,length=120),3)) |> 
  tidyr::drop_na() |> 
  as.matrix()

# static factors
d_static <- d |> 
  sf::st_drop_geometry() |> 
  dplyr::select(-c(id,day1.00:day5.23))

# building predictor list 
precipitation <- d_prec

# setting time steps
time <- (seq(1, 6, length = 120)) ## evaluation points
Time <- matrix(time, 3233, 120, byrow=TRUE) 

# fitting list
predlist <- list(Time=Time, 
                 precipitation=precipitation,
                 bin = as.numeric(as.vector(d_static$bin)),
                 doy = as.vector(d_static$doy),
                 month = as.vector(d_static$month),
                 year = as.vector(d_static$year),
                 elevation = as.vector(d_static$elevation),
                 slope = as.vector(d_static$slope),
                 aspect = as.vector(d_static$aspect),
                 standheight = as.vector(d_static$standheight),
                 normheight = as.vector(d_static$normheight),
                 concavity = as.vector(d_static$concavity),
                 convexity = as.vector(d_static$convexity),
                 curvature = as.vector(d_static$generalcurv),
                 convergence = as.vector(d_static$convergence),
                 twi = as.vector(d_static$twi),
                 tpi = as.vector(d_static$tpi),
                 swi = as.vector(d_static$swi),
                 catchment_area = as.vector(d_static$catchment_area),
                 mean_prec = as.vector(d_static$mean_prec),
                 landcover = as.vector(d_static$landcover),
                 geology = as.vector(d_static$geology),
                 forest = as.vector(d_static$forest),
                 geomorphon = as.vector(d_static$geomorphonre))


# FIT ---------------------------------------------------------------------

# cumulative precipitation 
fit <- refund::pfr(bin ~ af(precipitation, k=c(5,5)) +
                     s(slope, k=5) +
                     s(standheight, k=3) +
                     s(mean_prec, k=3) +
                     s(doy, bs="cc") +
                     as.factor(forest) +
                     as.factor(geology),
                   data=predlist, family="binomial");summary(fit)

# fitting performance (AUROC)
unlist(pROC::roc(predlist$bin, as.numeric(predict(fit, newdata=predlist, type="response", PredOutRange=T)), auc=T))$auc



# PLOTS -------------------------------------------------------------------

# interaction effect: precipitation and time
vis.gam(fit, view=c("precipitation.tmat","precipitation.omat"), plot.type="contour",
        too.far=0.2, n.grid=200, color="cm")

# partial effect plots for nonlinear terms
par(mfrow = c(2,2));par(pty="s")
plot(fit, select=4, trans=plogis, scale=-1, seWithMean=T, shade=T, shade.col="grey", ylab="", xlab="Mean annual precipitation")
plot(fit, select=5, trans=plogis, scale=-1, seWithMean=T, shade=T, shade.col="grey", ylab="", xlab="Day of year")
plot(fit, select=2, trans=plogis, scale=-1, seWithMean=T, shade=T, shade.col="grey", ylab="", xlab="Slope steepness")
plot(fit, select=3, trans=plogis, scale=-1, seWithMean=T, shade=T, shade.col="grey", ylab="", xlab="Standard height")

# arranging plots for forest and lithology factors
fe <- dplyr::tibble(name = c("No Forest", "Forest", "L1", "L2", "L3", "L4", "L5"),
                    group = c("land cover", "land cover", rep("lithology", 5)),
                    rc = c(0, summary(fit)$p.coeff[2], 0, summary(fit)$p.coeff[3:6]),
                    se = c(0, summary(fit)$se[2], 0, summary(fit)$se[3:6])) |> 
  dplyr::mutate_at(c('name', 'group'), as.factor)


# plots forest and lithology
gridExtra::grid.arrange(
dplyr::filter(fe, group == "land cover") |>
  dplyr::mutate(name = forcats::fct_relevel(name, 'No Forest', 'Forest')) |> 
  ggplot(aes(x = name, y = rc))+
  geom_errorbar(aes(ymin = rc-2*se, ymax = rc+2*se), linewidth=0.75, width = 0.5) +
  geom_point(aes(x = name, y = rc), size = 2, color = 'red') +
  theme_plot(),

dplyr::filter(fe, group == "lithology") |>
  ggplot(aes(x = name, y = rc))+
  geom_errorbar(aes(ymin = rc-2*se, ymax = rc+2*se), linewidth=0.75, width = 0.5) +
  geom_point(aes(x = name, y = rc), size = 2, color = 'red') +
  theme_plot(),
ncol = 2)


# VARIABLE IMPORTANCE -----------------------------------------------------

# setting formula
formula <- bin ~
  af(precipitation, k=c(5,5)) +
  s(slope, k=5) +
  s(standheight, k=3) +
  s(mean_prec, k=3) +
  s(doy, bs="cc") +
  as.factor(forest) +
  as.factor(geology)

#  setting baseline
mod <- refund::pfr(formula=formula, family=binomial, data=predlist);summary(mod)
formula <-mod$formula
formula.terms <- labels(terms(formula))
formula.vars <- all.vars(formula)[-3]
prop.dev <- rep(NA, length(formula.terms))

# setting null model
mod.null <- mgcv::gam(formula(paste0(formula.vars[1], '~1')), family=binomial, data=predlist);summary(mod.null)

# iterative to extract the deviance explained for every term
for(i in 1:length(formula.terms)){
  if (i == 1){
    formula.i <-  formula(paste0(formula.vars[1], "~", paste0(formula.terms[-i], collapse = "+")))
    mod.i <-  refund::pfr(formula.i, sp = mod$sp[-c(1,2)], family = binomial, data = predlist)
  } else {
    formula.i <-  formula(paste0(formula.vars[1], "~", paste0("af(precipitation, k=c(5,5)) + "), paste0(formula.terms[-c(1,i)], collapse = "+")))
    mod.i <-  refund::pfr(formula.i, sp = mod$sp[-c(i+1)], family = binomial, data = predlist)
  }
  prop.dev[i] <-  ((deviance(mod.i)-deviance(mod))/deviance(mod.null))
}

# arranging results for ploting
formula.vars <- formula.vars[-c(1, 2, 3)]
formula.vars <- c('precipitation', formula.vars)
dev.explain <- dplyr::tibble(name = formula.vars,
                             deviance = prop.dev) |> 
  dplyr::arrange(deviance) |> 
  dplyr::mutate(name = factor(name, levels=name)) |> 
  dplyr::mutate(n = 1:7)


# barplot of variable importance
ggplot(dev.explain, aes(x=name, y=deviance, color=name)) +
  geom_segment(aes(xend=name, y=0, yend=deviance), linewidth=1)+ geom_point(color="black") +
  scale_color_manual(values = (colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(7))) +
  ylab("Proportion of deviance explained") +
  xlab("") +  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.3), breaks = seq(0, 0.3, by = 0.05)) +
  coord_flip() + theme_plot()



