library(urca)
library(forecast)
library(tidyverse)
library(vars)
library(tsDyn)
library(ggplot2)
library(patchwork)
library(cowplot)
library(gridExtra)
library(pdftools)

#library(devtools)

## Function for extracting dataframe from object (irfs)

extract_varirf <- function(...){
  
  varirf_object <- list(...) #list one or more varirf input objects
  
  get_vec_length <- function(list_item){nrow(list_item[[1]][[1]])}
  
  if (!("varirf" %in% mapply(class, varirf_object))){
    stop("this function only accepts 'varirf' class objects")
  }
  
  if (length(unique(mapply(class, varirf_object)))!=1){
    stop("all input items must be 'varirf' class objects")
  }    
  if (length(unique(mapply(get_vec_length, varirf_object)))!=1){
    stop("all irf vectors must have the same length")   
  }  
  
  period <- as.data.frame(0:(nrow(varirf_object[[1]][[1]][[1]])-1)) 
  names(period) <- "period"
  
  for (l in 1:length(varirf_object)){
    for (i in 1:3){
      for (j in 1:dim(varirf_object[[l]][[i]][[1]])[2]){
        for (k in 1:length(varirf_object[[l]][[1]])){
          temp_colname <- paste(names(varirf_object[[l]][i]), #vector type (irf, lower, or upper)
                                names(varirf_object[[l]][[i]])[k], #impulse name
                                colnames(varirf_object[[l]][[i]][[k]])[j], #response name
                                sep = "_")
          
          temp <- as.data.frame(varirf_object[[l]][[i]][[k]][, j]) #extracts the vector
          
          names(temp) <- temp_colname #add the column name (vectortype_impulse_reponse)
          period <- cbind(period, temp) 
        }
        
      }
    }
  }
  names(period) <- tolower(names(period))
  return(period)
}

###2007M10 2019M12

setwd("C:/R_data")
dataset =read.csv("data_ret_mon_pol_R.csv")

data_post =subset(dataset, Period >= 123 & Period <= 268) #2007 nov - 2019 dec
data_ext =subset(dataset, Period >= 29 & Period <= 133) #1999 jan - 2020 dec

lev_post <- 100*ts(data_post$gamma_ln, frequency = 12)
lev_ext <- 100*ts(data_ext$lev_ln, frequency = 12)

hicp_post <- ts(data_post$hicp_ln, frequency = 12)
hicp_ext <- ts(data_ext$hicp_ln, frequency = 12)

r_post <- ts(data_post$real, frequency = 12)
r_ext <- ts(data_ext$real_3m, frequency = 12)

D_dc_ext <- ts(data_ext$D_dc, frequency = 12)
D_zlb_ext <- ts(data_ext$D_zlb, frequency = 12)

D_dc_post <- ts(data_post$D_dc, frequency = 12)
D_zlb_post <- ts(data_post$D_zlb, frequency = 12)

### Cholesky identification

dset_post <- cbind(hicp_post,r_post,lev_post) 
dset_ext <- cbind(hicp_ext,r_ext,lev_ext)

dset_dum_post <- cbind(D_dc_post,D_zlb_post)
dset_dum_ext <- cbind(D_dc_ext,D_zlb_ext)

lagselect_post <- VARselect(dset_post, lag.max = 12, type = "const")
lagselect_ext <- VARselect(dset_ext, lag.max = 12, type = "const")

lagselect_post$selection
lagselect_ext$selection

ctest1t_post <- ca.jo(dset_post, type = "trace", ecdet = "trend", K = 4, spec="longrun", season = 12, dumvar = dset_dum_post)
summary(ctest1t_post)
ctest2t_post <- ca.jo(dset_post, type = "trace", ecdet = "trend", K = 4, spec="longrun", season = 12, dumvar = dset_dum_post)
summary(ctest2t_post)

ctest1t_ext <- ca.jo(dset_ext, type = "eigen", ecdet = "trend", K = 4, spec="longrun", season = 12)
summary(ctest1t_ext)
ctest2t_ext <- ca.jo(dset_ext, type = "eigen", ecdet = "trend", K = 4, spec="longrun", season = 12)
summary(ctest2t_ext)

Model_post <- lineVar(
  dset_post,
  4,
  r = 2,
  include = c("none"),
  model = c("VECM"),
  #I = c("level"),
  beta = NULL,
  estim = c("ML"),
  LRinclude = c("both"),
  exogen = dset_dum_post
)

Model_ext <- lineVar(
  dset_ext,
  4,
  r = 2,
  include = c("none"),
  model = c("VECM"),
  #I = c("level"),
  beta = NULL,
  estim = c("ML"),
  LRinclude = c("both"),
)

summary(Model_post)
summary(Model_ext)

Model_VAR_post <- vec2var(ctest1t_post, r=2)
Model_VAR_ext <- vec2var(ctest1t_ext, r=2)

## Output figures

Irfs_lev_post <- irf(Model_VAR_post, impulse = "lev_post", n.ahead = 24, boot = TRUE)
Irfs_lev_ext <- irf(Model_VAR_ext, impulse = "lev_ext", n.ahead = 24, boot = TRUE)

Irfs_r_post <- irf(Model_VAR_post, impulse = "r_post", n.ahead = 24, boot = TRUE)
Irfs_r_ext <- irf(Model_VAR_ext, impulse = "r_ext", n.ahead = 24, boot = TRUE)

Irfs_hicp_post <- irf(Model_VAR_post, impulse = "hicp_post", n.ahead = 24, boot = TRUE)
Irfs_hicp_ext <- irf(Model_VAR_ext, impulse = "hicp_ext", n.ahead = 24, boot = TRUE)

## Leverage IRFs

single_varirf_lev_ext <- extract_varirf(Irfs_lev_ext)
single_varirf_lev_post <- extract_varirf(Irfs_lev_post)
head(single_varirf_lev_ext)
head(single_varirf_lev_post)

irf_lev_ext_r_ext <- single_varirf_lev_ext %>% 
  ggplot(aes(x=period, y=irf_lev_ext_r_ext , ymin=lower_lev_ext_r_ext , ymax=upper_lev_ext_r_ext)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  ggtitle("Response to Leverage")+
  ylab("")+
  xlab("") +
  #ylim(-0.1,+0.25)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_lev_ext_r_ext

irf_lev_post_r_post <- single_varirf_lev_post %>% 
  ggplot(aes(x=period, y=irf_lev_post_r_post , ymin=lower_lev_post_r_post , ymax=upper_lev_post_r_post)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  ggtitle("Response to Leverage")+
  ylab("")+
  xlab("") +
  #ylim(-0.1,+0.25)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_lev_post_r_post

irf_lev_ext_hicp_ext <- single_varirf_lev_ext %>% 
  ggplot(aes(x=period, y=irf_lev_ext_hicp_ext , ymin=lower_lev_ext_hicp_ext , ymax=upper_lev_ext_hicp_ext)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Leverage")+
  ylab("")+
  xlab("") +
  #ylim(-0.05,+0.075)+ 
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_lev_ext_hicp_ext

irf_lev_post_hicp_post <- single_varirf_lev_post %>% 
  ggplot(aes(x=period, y=irf_lev_post_hicp_post , ymin=lower_lev_post_hicp_post , ymax=upper_lev_post_hicp_post)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Leverage")+
  ylab("")+
  xlab("") +
  # #ylim(-0.05,+0.075)+
  # #ylim(-0.05,+0.3)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_lev_post_hicp_post

irf_lev_ext_lev_ext <- single_varirf_lev_ext %>% 
  ggplot(aes(x=period, y=irf_lev_ext_lev_ext, ymin=lower_lev_ext_lev_ext, ymax=upper_lev_ext_lev_ext )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Leverage")+
  ylab("")+
  xlab("") +
  #ylim(-0.15,+1.1)+ 
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_lev_ext_lev_ext

irf_lev_post_lev_post <- single_varirf_lev_post %>% 
  ggplot(aes(x=period, y=irf_lev_post_lev_post, ymin=lower_lev_post_lev_post, ymax=upper_lev_post_lev_post )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Leverage")+
  ylab("")+
  xlab("") +
  #ylim(-0.15,+1.1)+ 
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_lev_post_lev_post

irf_lev <- (irf_lev_ext_r_ext-irf_lev_post_r_post )/(irf_lev_ext_hicp_ext - irf_lev_post_hicp_post)/(irf_lev_ext_lev_ext-irf_lev_post_lev_post)

irf_lev

## HICP IRFs

single_varirf_hicp_ext <- extract_varirf(Irfs_hicp_ext)
head(single_varirf_hicp_ext)
single_varirf_hicp_post <- extract_varirf(Irfs_hicp_post)
head(single_varirf_hicp_post)


irf_hicp_ext_r_ext  <- single_varirf_hicp_ext %>% 
  ggplot(aes(x=period, y=irf_hicp_ext_r_ext  , ymin=lower_hicp_ext_r_ext   , ymax=upper_hicp_ext_r_ext )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  ggtitle("Response to Core HICP")+
  ylab("")+
  xlab("") +
  #ylim(-0.12,+0.17)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_hicp_ext_r_ext 

irf_hicp_post_r_post  <- single_varirf_hicp_post %>% 
  ggplot(aes(x=period, y=irf_hicp_post_r_post  , ymin=lower_hicp_post_r_post , ymax=upper_hicp_post_r_post)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  ggtitle("Response to Core HICP")+
  ylab("")+
  xlab("") +
  #ylim(-0.12,+0.17)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_hicp_post_r_post 

irf_hicp_ext_hicp_ext <- single_varirf_hicp_ext %>% 
  ggplot(aes(x=period, y=irf_hicp_ext_hicp_ext  , ymin=lower_hicp_ext_hicp_ext , ymax=upper_hicp_ext_hicp_ext )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Core HICP")+
  ylab("")+
  xlab("") +
  #ylim(-0.025,+0.11)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_hicp_ext_hicp_ext

irf_hicp_post_hicp_post <- single_varirf_hicp_post %>% 
  ggplot(aes(x=period, y=irf_hicp_post_hicp_post  , ymin=lower_hicp_post_hicp_post , ymax=upper_hicp_post_hicp_post)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Core HICP")+
  ylab("")+
  xlab("") +
  #ylim(-0.025,+0.11)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_hicp_post_hicp_post

irf_hicp_ext_lev_ext <- single_varirf_hicp_ext %>% 
  ggplot(aes(x=period, y=irf_hicp_ext_lev_ext, ymin=lower_hicp_ext_lev_ext , ymax=upper_hicp_ext_lev_ext   )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Core HICP")+
  ylab("")+
  xlab("") +
  #ylim(-0.28,+0.5)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_hicp_ext_lev_ext 

irf_hicp_post_lev_post <- single_varirf_hicp_post %>% 
  ggplot(aes(x=period, y=irf_hicp_post_lev_post, ymin=lower_hicp_post_lev_post , ymax=upper_hicp_post_lev_post )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Core HICP")+
  ylab("")+
  xlab("") +
  #ylim(-0.28,+0.5)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_hicp_post_lev_post 

irf_hicp <- (irf_hicp_ext_r_ext-irf_hicp_post_r_post )/(irf_hicp_ext_hicp_ext - irf_hicp_post_hicp_post)/(irf_hicp_ext_lev_ext-irf_hicp_post_lev_post)

irf_hicp

## Real int IRFs

single_varirf_r <- extract_varirf(Irfs_r_ext, Irfs_r_post)
head(single_varirf_r)

irf_r_ext_r_ext  <- single_varirf_r %>% 
  ggplot(aes(x=period, y=irf_r_ext_r_ext  , ymin=lower_r_ext_r_ext   , ymax=upper_r_ext_r_ext )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ylim(0,+0.38)+ 
  scale_x_continuous(breaks=seq(0, 24, 6))+
  ggtitle("Response to Real Interest")+
  ylab("Real Interest Rate")+
  xlab("") +
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_r_ext_r_ext 

irf_r_post_r_post  <- single_varirf_r %>% 
  ggplot(aes(x=period, y=irf_r_post_r_post  , ymin=lower_r_post_r_post   , ymax=upper_r_post_r_post )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ylim(0,+0.38)+ 
  ggtitle("Response to Real Interest")+
  ylab("Real Interest Rate")+
  xlab("") +
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_r_post_r_post 

irf_r_ext_hicp_ext <- single_varirf_r %>% 
  ggplot(aes(x=period, y=irf_r_ext_hicp_ext  , ymin=lower_r_ext_hicp_ext , ymax=upper_r_ext_hicp_ext )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Real Interest")+
  ylab("Core HICP")+
  xlab("") +
  #ylim(-0.085,+0.03)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_r_ext_hicp_ext

irf_r_post_hicp_post <- single_varirf_r %>% 
  ggplot(aes(x=period, y=irf_r_post_hicp_post  , ymin=lower_r_post_hicp_post , ymax=upper_r_post_hicp_post )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  # ##ylim(-0.45,+0.45)+ 
  #ggtitle("Response to Real Interest")+
  ylab("Core HICP")+
  xlab("") +
  #ylim(-0.085,+0.03)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_r_post_hicp_post

irf_r_ext_lev_ext <- single_varirf_r %>% 
  ggplot(aes(x=period, y=irf_r_ext_lev_ext, ymin=lower_r_ext_lev_ext , ymax=upper_r_ext_lev_ext   )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Real Interest")+
  ylab("Leverage Ratio")+
  xlab("") +
  #ylim(-0.6,+0.55)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

#irf_r_ext_lev_ext 


irf_r_post_lev_post <- single_varirf_r %>% 
  ggplot(aes(x=period, y=irf_r_post_lev_post, ymin=lower_r_post_lev_post , ymax=upper_r_post_lev_post   )) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="light blue", alpha=0.2) +
  geom_line() +
  theme_light() +
  #ggtitle("Response to Real Interest")+
  ylab("Leverage Ratio")+
  xlab("") +
  #ylim(-0.6,+0.55)+
  scale_x_continuous(breaks=seq(0, 24, 6))+
  theme(plot.title = element_text(size = 10, hjust=0.5),
        axis.title.y = element_text(size=10))

irf_post <- (irf_r_post_r_post | irf_lev_post_r_post | irf_hicp_post_r_post) / 
  (irf_r_post_hicp_post | irf_lev_post_hicp_post | irf_hicp_post_hicp_post) /
  (irf_r_post_lev_post | irf_lev_post_lev_post | irf_hicp_post_lev_post)

irf_post

pdf(file="irf_post.pdf")
irf_post
dev.off()

irf_ext <- (irf_r_ext_r_ext | irf_lev_ext_r_ext | irf_hicp_ext_r_ext) / 
  (irf_r_ext_hicp_ext | irf_lev_ext_hicp_ext | irf_hicp_ext_hicp_ext) /
  (irf_r_ext_lev_ext | irf_lev_ext_lev_ext | irf_hicp_ext_lev_ext)

irf_ext

pdf(file="irf_ext.pdf")
irf_ext
dev.off()