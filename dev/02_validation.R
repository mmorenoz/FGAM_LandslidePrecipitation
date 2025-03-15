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
list.packages <- c("mgcv", "refund" ,"sf", "tidyverse", "sperrorest")
vapply(list.packages, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)
remove(list.packages)

# ggplot theme
theme_plot <- function(){
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.box.background = element_rect(color="black", linewidth=0.3),
        legend.box.margin = margin(1, 1, 1, 1),
        panel.grid.minor = element_blank(),
        axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14, vjust= 2),
        axis.title.x  = element_text(size=14, vjust=-2),
        legend.position = c(0.11,0.9),
        plot.margin = unit(c(1,2,1,1), "lines"),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5))
}


# LOADING DATA ------------------------------------------------------------

# load cum precipitation
d <- readRDS("./dat/d_5_0cumforward_5day_t0123.Rds") |> 
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
                 slope = as.vector(d_static$slope),
                 standheight = as.vector(d_static$standheight),
                 mean_prec = as.vector(d_static$mean_prec),
                 geology = as.factor(d_static$geology),
                 forest = as.factor(d_static$forest))


# VALIDATION --------------------------------------------------------------

# separate X and Y coordinates for plot
d <- d |>  
  dplyr::mutate(x = sf::st_coordinates(d)[,1]) |>  
  dplyr::mutate(y = sf::st_coordinates(d)[,2]) |>  
  st_drop_geometry()

# formula
formula.cv <- bin ~
  af(precipitation, k=c(5,5)) +
  s(slope, k=5) +
  s(standheight, k=3) +
  s(mean_prec, k=3) +
  s(doy, bs="cc") +
  as.factor(forest) +
  as.factor(geology)
  

## random cv ---------------------------------------------------------------

# setting up loop for random cross-validation
fold <- 10
repetition <- 10

# create random cv partition
set.seed(1)
partition <- sperrorest::partition_cv(d, nfold = fold, repetition = repetition, seed1 = 123) 
myroc_cv <- lapply(my.list<-vector(mode = 'list', 10),function(x) x<-vector(mode='list', 10))

# loop for validation
for (i in 1:repetition){
  id.hold <- partition[[i]]
  for (j in 1:fold){
    id.holdout <- id.hold[[j]]$test
    predlist.test <- lapply(predlist, `[`, id.holdout)
    predlist.test$precipitation <- d_prec[id.holdout,]
    predlist.train <- lapply(predlist, `[`, -id.holdout) 
    predlist.train$precipitation <- d_prec[-id.holdout,]
    fit <- refund::pfr(formula.cv, data=predlist.train, family="binomial")
    predlist.test$prediction <- predict(fit, newdata=predlist.test, type="response")
    myroc <- unlist(pROC::roc(response=predlist.test$bin, predictor=predlist.test$prediction, auc=T))$auc
    myroc_cv[[i]][[j]] <- myroc
  }
}

# clean environment
rm(predlist.train, predlist.test, id.hold, id.holdout, bin, slope,
   standheight, mean_prec, forest, geology, doy, myroc, i, j, Time)

# store performances in data frame
myroc_cv = dplyr::as_tibble(do.call(rbind, lapply(myroc_cv, unlist))) |> 
  tidyr::gather(key="repetition", value="auc") |> 
  dplyr::mutate(repetition = rep(1:10, each=10))


## spatial cv ---------------------------------------------------------------

# setting up loop for spatial cross-validation
fold <- 10
repetition <- 10

# create spatial scv partition
set.seed(1)
partition.s <- sperrorest::partition_kmeans(d, nfold = fold, repetition = repetition, seed1 = 123) 
myroc_scv <- lapply(my.list<-vector(mode = 'list',10),function(x) x<-vector(mode='list',10))

# loop for validation
for (i in 1:repetition){
  id.hold <- partition.s[[i]]
  for (j in 1:fold){
    id.holdout <- id.hold[[j]]$test
    predlist.test <- lapply(predlist, `[`, id.holdout)
    predlist.test$precipitation <- d_prec[id.holdout,]
    predlist.train <- lapply(predlist, `[`, -id.holdout) 
    predlist.train$precipitation <- d_prec[-id.holdout,]
    fit <- refund::pfr(formula.cv, data=predlist.train, family="binomial")
    predlist.test$prediction <- predict(fit, newdata=predlist.test, type="response")
    myroc <- unlist(pROC::roc(response=predlist.test$bin, predictor=predlist.test$prediction, auc=T))$auc
    myroc_scv[[i]][[j]] <- myroc
  }
}

# clean environment
rm(predlist.train, predlist.test, id.hold, id.holdout, bin, slope,
   standheight, mean_prec, forest, geology, doy, myroc, i, j, Time)

# store performances in data frame
myroc_scv <- dplyr::as_tibble(do.call(rbind, lapply(myroc_scv, unlist))) |>  
  tidyr::gather(key="repetition", value="auc") |> 
  dplyr::mutate(repetition = rep(1:10, each=10))


## factor cv ---------------------------------------------------------------

### lithology --------------------------------------------------------------

summary(predlist$geology)

# create partitions
set.seed(1)
partition.geo <- sperrorest::partition_factor(d, fac = "geology") 

# setting up loop for cross-validation
fold <- length(partition.geo[[1]])
myroc_gcv <- c()

# loop for validation
for (i in 1:fold){
  id.holdout <- partition.geo[[1]][[i]]$test 
  predlist.test <- lapply(predlist, `[`, id.holdout)
  predlist.test$precipitation <- d_prec[id.holdout,]
  predlist.train <- lapply(predlist, `[`, -id.holdout) 
  predlist.train$precipitation <- d_prec[-id.holdout,]
  fit <- refund::pfr(formula.cv, data=predlist.train, family="binomial", drop.unused.levels=F)
  predlist.test$prediction <- predict(fit, newdata=predlist.test, type="response", newdata.guaranteed=TRUE)
  myroc <- unlist(pROC::roc(response=predlist.test$bin, predictor=predlist.test$prediction, auc=T))$auc
  myroc_gcv[i] <- myroc
}

# clean environment
rm(predlist.train, predlist.test, id.holdout,bin, slope,
   doy, fold, forest, geology, myroc, mean_prec, i, Time)

# store performances in tibble
myroc_gcv <- dplyr::tibble(fold = c(1:5),
                           geology = c(
                             "Crystalline",
                             "Porphyry",
                             "Sedimentary",
                             "Plutonic",
                             "Calcschist"), 
                           auc = myroc_gcv) |> 
  dplyr::mutate(auc = as.numeric(auc)) |> 
  dplyr::mutate(fold = as.integer(fold)) |> 
  dplyr::arrange(fold)


### month --------------------------------------------------------------

# convert to factor
d <- dplyr::mutate(d, month = as.factor(month))
summary(d$month)

# create partition
set.seed(1)
partition.month <- sperrorest::partition_factor(d, fac = "month") 
plot(partition.month, d, cex = 0.01, pch = 19)

# setting up loop for cross-validation
fold <- length(partition.month[[1]])
myroc_mcv <- c()

# loop for validation
for (i in 1:fold){
  id.holdout <- partition.month[[1]][[i]]$test
  predlist.test <- lapply(predlist, `[`, id.holdout)
  predlist.test$precipitation <- d_prec[id.holdout,]
  predlist.train <- lapply(predlist, `[`, -id.holdout) 
  predlist.train$precipitation <- d_prec[-id.holdout,]
  fit <- refund::pfr(formula.cv, data=predlist.train, family="binomial", drop.unused.levels=F)
  predlist.test$prediction <- predict(fit, newdata=predlist.test, type="response", newdata.guaranteed=TRUE)
  myroc <- unlist(pROC::roc(response=predlist.test$bin, predictor=predlist.test$prediction, auc=T))$auc
  myroc_mcv[i] <- myroc
}

# clean environment
rm(predlist.train, predlist.test, id.holdout,bin, slope,
   doy, fold, forest, geology, myroc, mean_prec, i, Time)

# store performances in data frame
myroc_mcv <- dplyr::tibble(fold = seq(1,12,1),
                           month = seq(1,12,1), 
                           auc = myroc_mcv)|> 
  dplyr::mutate(auc = as.numeric(auc)) |> 
  dplyr::mutate(fold = as.integer(fold)) |> 
  dplyr::arrange(month)


### year --------------------------------------------------------------

# convert to factor
d <- dplyr::mutate(d, year = as.factor(year))
summary(d$year)

# create partition
set.seed(1)
partition.year <- sperrorest::partition_factor(d, fac = "year") 
plot(partition.year, d, cex = 0.01, pch = 19)

# setting up loop for cross-validation
fold <- length(partition.year[[1]])
myroc_ycv <- c()

# loop for validation
for (i in 1:fold){
  id.holdout <- partition.year[[1]][[i]]$test
  predlist.test <- lapply(predlist, `[`, id.holdout)
  predlist.test$precipitation <- d_prec[id.holdout,]
  predlist.train <- lapply(predlist, `[`, -id.holdout) 
  predlist.train$precipitation <- d_prec[-id.holdout,]
  fit <- refund::pfr(formula.cv, data=predlist.train, family="binomial", drop.unused.levels=F)
  predlist.test$prediction <- predict(fit, newdata=predlist.test, type="response", newdata.guaranteed=TRUE)
  myroc <- unlist(pROC::roc(response=predlist.test$bin, predictor=predlist.test$prediction, auc=T))$auc
  myroc_ycv[i] <- myroc
}

# clean environment
rm(predlist.train, predlist.test, id.holdout,bin, slope,
   doy, fold, forest, geology, myroc, mean_prec, i, Time)

# store performances in data frame
myroc_ycv = dplyr::tibble(fold = seq(1,10,1),
                          year = c(2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021), 
                          auc = myroc_ycv) |> 
  dplyr::mutate(auc = as.numeric(auc)) |> 
  dplyr::mutate(fold = as.integer(fold)) |> 
  dplyr::arrange(year)


# store environment data
save.image("./dat/CV.RData")


# PLOTS -------------------------------------------------------------------

# loading data
load("./dat/CV.RData")

# stats
#### RANDOM AND SPATIAL ----
pdf("./plt/03_rcdandscv.pdf", width = 11, height = 8, paper="a4r")
par(mfrow=c(1,2))
boxplot(myroc_cv$auc,
        ylim=c(0.8, 1),
        boxwex=0.2,
        main = paste0(paste0("IQR= "), paste0(round(IQR(myroc_cv$auc),4)),
                      paste0(" median AUC= "), paste0(round(median(myroc_cv$auc), 3))))
boxplot(myroc_scv$auc,
        ylim=c(0.8, 1),
        boxwex=0.2,
        main = paste0(paste0("IQR= "), paste0(round(IQR(myroc_scv$auc),4)),
                      paste0(" median AUC= "), paste0(round(median(myroc_scv$auc), 3))))
dev.off()
par(mfrow=c(1,1))


## GEOLOGY ----

# loop to extract the number of lansaldies per fold
geology_frame <- dplyr::as_tibble(cbind(fold = rep(NA, 5), test = rep(NA, 5), landslides = rep(NA, 5)))

for (i in 1:5){
  id.holdout <- partition.geo[[1]][[i]]$test 
  test <- lapply(predlist, `[`, id.holdout)$bin
  geology_frame$fold[i] <- i
  geology_frame$test[i] <- length(test)
  geology_frame$landslides[i] <- length(test[test == 1])
  
}

# join number of landslides with auc dataframe
myroc_gcv <- myroc_gcv |>  
  dplyr::left_join(geology_frame, by="fold") |>  
  dplyr::mutate(geology = as.factor(geology)) |>  
  dplyr::mutate(geology = ordered(geology, levels = c("Crystalline", "Porphyry", "Sedimentary", "Plutonic", "Calcschist")))

# plot
ggplot(myroc_gcv, aes(x=geology, y=auc)) + geom_bar(stat="identity", width=0.75) +
  xlab("Lithology") + ylab("AUC") + ggtitle(paste("mean AUC:", round(mean(myroc_gcv$auc),3))) +
  coord_cartesian(ylim=c(0.60,1)) +
  geom_text(aes(label=test), vjust=4, color="white", size=3.5)+
  geom_text(aes(label=landslides), vjust=1.6, color="red", size=3.5)+
  theme_plot()
ggsave("./plt/04_geocv.pdf", width = 7, height = 8, dpi = 300, units = "in", device='pdf')


## MONTH ----
month_frame <- dplyr::as_tibble(cbind(fold = rep(NA, 12), test = rep(NA, 12), landslides = rep(NA, 12)))

# loop to extract the number of landslides per fold
for (i in 1:12){
  id.holdout <- partition.month[[1]][[i]]$test 
  test <- lapply(predlist, `[`, id.holdout)$bin
  month_frame$fold[i] <- i
  month_frame$test[i] <- length(test)
  month_frame$landslides[i] <- length(test[test == 1])
}

# join number of landslides with auc dataframe
myroc_mcv <- myroc_mcv |>  
  dplyr::left_join(month_frame, by="fold") |>  
  dplyr::mutate(month = as.factor(month))

# plot
ggplot(myroc_mcv, aes(x=as.factor(month), y=auc)) + geom_bar(stat="identity", width=0.75) +
  xlab("Month") + ylab("AUC") + ggtitle(paste("mean AUC:", round(mean(myroc_mcv$auc),3))) +
  coord_cartesian(ylim=c(0.60,1)) +
  geom_text(aes(label=test), vjust=4, color="white", size=3.5)+
  geom_text(aes(label=landslides), vjust=1.6, color="red", size=3.5)+
  theme_plot()
ggsave("./plt/05_month.pdf", width = 7, height = 8, dpi = 300, units = "in", device='pdf')

## YEAR ----
year_frame <- as.data.frame(cbind(fold = rep(NA, 10), test = rep(NA, 10), landslides = rep(NA, 10)))

# loop to extract the number of landslides per fold
for (i in 1:10){
  id.holdout <- partition.year[[1]][[i]]$test 
  test <- lapply(predlist, `[`, id.holdout)$bin
  year_frame$fold[i] <- i
  year_frame$test[i] <- length(test)
  year_frame$landslides[i] <- length(test[test == 1])
}

# join number of landslides with auc dataframe
myroc_ycv <- myroc_ycv |> 
  dplyr::left_join(year_frame, by="fold") |> 
  dplyr::mutate(year = as.factor(year))

# plot
ggplot(myroc_ycv, aes(x=as.factor(year), y=auc)) + geom_bar(stat="identity", width=0.75) +
  xlab("Year") + ylab("AUC") + ggtitle(paste("mean AUC:", round(mean(myroc_ycv$auc),3))) +
  coord_cartesian(ylim=c(0.60,1)) +
  geom_text(aes(label=test), vjust=4, color="white", size=3.5)+
  geom_text(aes(label=landslides), vjust=1.6, color="red", size=3.5)+
  theme_plot()
ggsave("./plt/06_yearcv.pdf", width = 7, height = 8, dpi = 300, units = "in", device='pdf')


