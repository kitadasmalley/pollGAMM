LD_LIBRARY_PATH=/home/stats/wolfch/lib:/cm/local/apps/gcc/5.1.0/lib64:/cm/shared/apps/gdal/2.1.1/lib:/cm/shared/apps/gdal/2.1.1/bin:/home/stats/wolfch/my_libs:/home/stats/wolfch/R_libs/3.5.1:/cm/shared/apps/R/3.5.1/lib64/R:/home/stats/wolfch/wood/curl/wolf_curl/lib:/cm/shared/apps/R/3.2.2/lib64/R:/cm/shared/apps/slurm/14.11.11/lib64/slurm:/cm/shared/apps/slurm/14.11.11/lib64:/cm/shared/apps/sge/2011.11p1/lib/linux-x64:/cm/local/apps/gcc/5.1.0/lib:/cm/local/apps/gcc/5.1.0/lib64
export LD_RUN_PATH=/cm/shared/apps/gdal/2.1.1/bin
export PKG_CONFIG_PATH=/home/stats/wolfch/lib/pkgconfig
export R_LIBS=/home/stats/wolfch/R_libs/3.5.1

/cm/shared/apps/R/3.5.1/bin/R

require(sf)
require(maps)
require(unixtools)
require(rstan)
library(tidyverse)
library(scales)
library(pscl)
library(forcats)
require(dslabs)
require(lubridate)
library(RCurl)
require(reshape2)
require(boot)
require(brms)
require(mgcv)
require(itsadug)
require(ggplot2)
require(RColorBrewer)
require(zoo)
require(plyr)
require(ggstance)
require(Hmisc)
require(gridExtra)
require(multcomp)
require(broom)
require(mgcv.helper) # devtools::install_github("samclifford/mgcv.helper")

options(mc.cores = 4)

setwd("/home/stats/wolfch/misc/Heather")
set.tempdir("/home/stats/wolfch/temp2")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

####################### 0. clean up the election data

df = read.csv("2016 Presidential Election Polls.csv", as.is=TRUE)
df$startdate = as.Date(df$startdate, "%m/%d/%Y")
df$enddate = as.Date(format(as.Date(df$enddate, "%m/%d/%Y"), "20%y/%m/%d"))
df$time = as.numeric(df$enddate)

df$Mode[df$Mode == "Telepone"] = "Telephone"
df = df[df$Mode %in% c("Telephone","Online"),]
df$Mode = factor(df$Mode)
df$pollster = factor(df$pollster, levels=c(unique(df$pollster),"new"))

df$state = gsub(" CD-.*","",df$state)
df$state.abb = state.abb[match(df$state, state.name)]
df$state.abb[df$state == "District of Columbia"] = "DC"
df$state = factor(df$state)

# add some state-level covariates
# df = data.frame(df, state.x77[match(df$state,rownames(state.x77)),]) # old version w/ built-in data
state_data = read.csv("ACS_16_dems.csv", as.is=TRUE)
df = data.frame(df, state_data[match(df$state,state_data$Geography),])

df$wts = 1/table(df$Mode)[df$Mode] # weighting each mode equally might be more important than weighting by pollster...
df$wts = df$wts/mean(df$wts) # since we include pollster as a random effect

# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/42MVDX
results = read.csv("1976-2016-president.csv", as.is=TRUE)
results = results[results$year == 2016 & results$candidate %in% c("Clinton, Hillary","Trump, Donald J.")
                  & ! results$writein & results$party %in% c("democrat","republican"),]
results = spread(results[,c("state","candidate","candidatevotes")],"candidate","candidatevotes")
trump_states = results$state[results$`Clinton, Hillary` < results$`Trump, Donald J.`]

####################### 0.5. additive mode effect models (easier to estimate parameters)

df[df$Mode == "Telephone", c("median_age","median_income","pct_below_poverty","pct_bachorhigher",
                          "sex_ratio_mw","pct_foreignborn","pct_white")] = 0
df$is_web = df$Mode == "Online"

# bam is super fast.  Stepwise variable selection based on p-values..
mod = bam(adjpoll_trump ~ s(time,by=state) + is_web + median_age + median_income + pct_below_poverty +
          pct_bachorhigher + sex_ratio_mw + pct_foreignborn + pct_white, data=df[df$state != "U.S.",],
          discrete=TRUE, nthreads=4)
summary(mod)
CI = data.frame(mgcv.helper:::confint.gam(mod, parm=c("is_webTRUE","median_age","median_income",
        "pct_below_poverty","pct_bachorhigher","sex_ratio_mw","pct_foreignborn","pct_white")))[,5:6]
apply(CI,1,function(x) paste0("(",paste(round(x,6),collapse=", "),")"))

mod = bam(adjpoll_trump ~ s(time,by=state) + is_web + median_age + median_income + pct_below_poverty +
          pct_foreignborn + pct_white, data=df[df$state != "U.S.",],
          discrete=TRUE, nthreads=4)
summary(mod)

mod = bam(adjpoll_trump ~ s(time,by=state) + is_web + median_age + median_income +
          pct_foreignborn + pct_white, data=df[df$state != "U.S.",],
          discrete=TRUE, nthreads=4)
summary(mod)

mod = bam(adjpoll_trump ~ s(time,by=state) + is_web + median_age + median_income +
          pct_white, data=df[df$state != "U.S.",],
          discrete=TRUE, nthreads=4)
summary(mod)

mod = bam(adjpoll_trump ~ s(time,by=state) + is_web + median_income + pct_white, data=df[df$state != "U.S.",],
          discrete=TRUE, nthreads=4)
summary(mod) # final model
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)    4.233e+01  1.809e-01 233.982  < 2e-16 ***
# is_webTRUE     1.009e+01  1.807e+00   5.584 2.60e-08 ***
# median_income -8.221e-05  2.093e-05  -3.928 8.80e-05 ***
# pct_white     -7.488e-02  1.773e-02  -4.223 2.49e-05 ***
data.frame(mgcv.helper:::confint.gam(mod, parm=c("is_webTRUE","median_income","pct_white")))

model.matrix(mod)[which(df[df$state != "U.S.",]$Mode == "Online")[1],] # see where mode effect columns are located

# https://stat.ethz.ch/pipermail/r-help/2012-August/320895.html
test_mat = 0 * model.matrix(mod)[1:51,]
test_mat[,2:4] = cbind(1, as.matrix(state_data[! state_data$Geography %in% c("United States","Puerto Rico"),
                                    c("median_income","pct_white")]))

# simultaneous CIs: https://www.rdocumentation.org/packages/multcomp/versions/1.4-10/topics/glht-methods
confint(glht(mod,linfct=test_mat)) # could compare to our other method (estimating prop. time O > T)

conf_df = tidy(confint(glht(mod,linfct=test_mat)))
conf_df$state = state_data$Geography[! state_data$Geography %in% c("United States","Puerto Rico")]
conf_df$state = factor(conf_df$state, levels=conf_df$state[order(conf_df$estimate)])
conf_df$trump = conf_df$state %in% trump_states
write.csv(conf_df, "output/conf_df.csv", row.names=FALSE)

p = ggplot(conf_df, aes(x=estimate,y=state,color=trump)) +
        geom_point() +
        geom_errorbarh(aes(xmin=conf.low, xmax=conf.high)) +
        theme_bw() +
  xlab("Mode effect (Online - Telephone)") +
  scale_color_manual(values=brewer.pal(3,"Set1")[2:1], guide=FALSE) +
  theme(axis.text=element_text(color="black"),
        axis.title.y=element_blank(),
        axis.ticks=element_line(color="black"),
        panel.border=element_rect(color="black")) +
  geom_vline(xintercept=0,linetype="dashed")

png("output/additive_mode_v2.png", width=7, height=7, units="in", res=300)
p
dev.off()

pdf("output/additive_mode_v2.pdf", width=7, height=7)
p
dev.off()

conf_df$not_sig = conf_df$conf.low < 0 & conf_df$conf.high > 0 # significant for 25 states; not sig. for 25 + DC

states = sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
conf_df$ID = tolower(conf_df$state)
states_merged = merge(states, conf_df)

p = ggplot(states_merged[! states_merged$not_sig,], aes(fill=estimate)) +
        geom_sf(data=states, fill="lightgray") +
        geom_sf() +
        scale_fill_gradient2(midpoint=0,
                             low = muted("purple"), mid = "white", high = muted("orange"),
                             guide=guide_colorbar(title="Mode effect (Online - Telephone)")) +
        theme_bw() +
        theme(axis.text=element_blank(),
              axis.ticks=element_blank(),
              legend.position="bottom")

png("output/state_map_v2.png", width=6, height=4, units="in", res=300)
p
dev.off()

pdf("output/state_map_v2.pdf", width=6, height=4)
p
dev.off()

####################### 1. test for model effect w/ compareML

for(response in c("adjpoll_trump","adjpoll_clinton")){ # national models (we can handle state models separately later)
        print(response)
        df$y = df[,response]

        national = gam(y ~ s(time,k=5) + s(pollster,bs="re"), weights=wts,
                       data=df[df$state == "U.S.",], method="ML")
        national_mode = gam(y ~ s(time,k=5,by=Mode) + s(pollster,bs="re"),  data=df[df$state == "U.S.",], method="ML")
        compareML(national, national_mode)
}

####################### 2. predict and plot (national, with and without mode)

pred = do.call("rbind", lapply(c("adjpoll_trump","adjpoll_clinton"), function(response){
        df$y = df[,response]

        national = gam(y ~ s(time) + s(pollster,bs="re"), data=df[df$state == "U.S.",], weights=wts,
                       method="ML",drop.unused.levels=FALSE)
        pred = data.frame(time = seq(min(df$time),max(df$time),length=1e2), Mode="Overall", pollster="new")
        pred$y = predict(national,newdata=pred,type="response")

        national_mode = gam(y ~ s(time,by=Mode) + s(pollster,bs="re"), data=df[df$state == "U.S.",],
                            method="ML",drop.unused.levels=FALSE)
        pred_mode = data.frame(expand.grid(time = seq(min(df$time),max(df$time),length=1e2), Mode=unique(df$Mode)),
                                         pollster="new")
        pred_mode$y = predict(national_mode,newdata=pred_mode,type="response")
        
        out = rbind(pred,pred_mode)
        out$variable = response
        return(out)
}))

poll_df = melt(df, measure.vars=c("adjpoll_trump","adjpoll_clinton"))

pred$Mode = factor(pred$Mode, levels=c("Overall","Online","Telephone"))
poll_df$Mode = factor(poll_df$Mode, levels=levels(pred$Mode))
pred$variable = factor(pred$variable, levels=c("adjpoll_clinton","adjpoll_trump"), labels=c("Clinton","Trump"))
poll_df$variable = factor(poll_df$variable, levels=c("adjpoll_clinton","adjpoll_trump"), labels=c("Clinton","Trump"))

p = ggplot(pred, aes(x=as.Date(time), y=value, color=Mode, size=Mode)) + 
        facet_wrap(~ variable) +
        theme_bw() +
        geom_point(data=poll_df[poll_df$state == "U.S.",], alpha=0.25,size=1,show.legend=FALSE) +
        geom_line(data=pred, aes(y=y,size=Mode)) +
  xlab("Poll date") +
  ylab("Polling percentage") +
  scale_color_manual(values=c("black",gg_color_hue(2)),drop=FALSE) +
  scale_size_manual(values=c(2,1,1),drop=FALSE) +
  theme(axis.text=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        panel.border=element_rect(color="black")) +
  geom_vline(xintercept=as.Date("7/28/2016","%m/%d/%Y"), linetype="dashed")

png("output/national_plot.png", width=7, height=4, units="in", res=300)
p
dev.off()

pdf("output/national_plot.pdf", width=7, height=4)
p
dev.off()

####################### 3. state-level model effect estimates (proportion time tel > web); link w/ covariates

model_mat = expand.grid(response=c("adjpoll_trump","adjpoll_clinton"),
                        state=names(which(table(df$state.abb,df$Mode)[,2] >= 5)), # at least 5 telephone polls
                        stringsAsFactors=FALSE)

df$y = df$adjpoll_trump
state_mod = brm(y ~ s(time,by=Mode), # + (1 | pollster),
          data = df[df$state.abb %in% "OR",], family = gaussian(),
          cores = 4, seed = 17, iter = 750, warmup = 500, thin = 1,
          silent=FALSE, refresh = 10, open_progress=TRUE) # fit model once and call with different data

state_results = do.call("rbind", lapply(1:nrow(model_mat), function(i){
        print(i)                                          
        newdata = df[df$state.abb %in% model_mat$state[i],]
        newdata$y = newdata[, model_mat$response[i]]
        state_mod = update(state_mod, newdata=newdata)
        #
        ms = marginal_smooths(state_mod,spaghetti=TRUE)
        ms_lines_long = attr(ms[[1]],"spaghetti")
        ms_lines_long$estimate__ = ms_lines_long$estimate__ + posterior_samples(state_mod,"b_Intercept")[ms_lines_long$sample__,1]
        ms_lines = spread(ms_lines_long, key="Mode", value="estimate__")
        ms_lines$estimate__ = ms_lines$Online - ms_lines$Telephone
        # for each sample, compute proportion of time that Online > Telephone
        prop_gt = tapply(ms_lines$Online > ms_lines$Telephone,ms_lines$sample__,mean)
        df_gt = data.frame(t(quantile(prop_gt, c(0.025,0.5,0.975))), Mode=NA)
        # also store predictions across time (for plotting); individuals modes version in comments
        #         pred = aggregate(estimate__ ~ time + Mode, data=ms_lines_long,
        #                 FUN=function(x) quantile(x,c(0.025,0.5,0.975)))
        pred = aggregate(estimate__ ~ time, data=ms_lines,
                         FUN=function(x) quantile(x,c(0.025,0.5,0.975)))
        pred = data.frame(pred, pred$estimate__)
        pred$estimate__ = NULL
        #
        out=data.frame(rbind.fill(pred,df_gt),
                response=model_mat$response[i],
                state=model_mat$state[i])
        return(out)
}))

state_results$response = factor(state_results$response,
        levels=c("adjpoll_trump","adjpoll_clinton"), labels=c("Trump","Clinton"))

saveRDS(state_results, "output/state_results.RDS")

####################### 4. plot/model state-level results

p = ggplot(state_results, aes(x=as.Date(time), y=X50., group=response, color=response)) +
        facet_wrap(~ state) +
        theme_bw() +
        geom_ribbon(aes(ymin=X2.5.,ymax=X97.5.,fill=response), color=NA, alpha=0.25) +
        geom_line() +
        coord_cartesian(ylim=c(-10,10)) +        
  xlab("Poll date") +
  ylab("Mode effect (Online - Telephone)") +
  scale_color_manual(values=brewer.pal(3,"Set1")[-3],
                     guide=guide_legend(title=NULL)) +
  scale_fill_manual(values=brewer.pal(3,"Set1")[-3],
                    guide=guide_legend(title=NULL)) +  
  theme(axis.text=element_text(color="black"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        panel.border=element_rect(color="black"),
        legend.position=1:0,
        legend.justification=1:0,
        legend.direction="horizontal") +
  geom_hline(yintercept=0, linetype="dashed")

png("output/state_lines_plot.png", width=7, height=7, units="in", res=300)
p
dev.off()

prop_gt = state_results[is.na(state_results$time),]
prop_gt$state = factor(prop_gt$state,
                       names(sort(tapply(prop_gt$X50.,prop_gt$state,mean))))

p = ggplot(prop_gt, aes(y=state, x=X50., color=response)) +
        geom_point(position=ggstance::position_dodgev(height=0.5)) +
        geom_errorbarh(aes(xmin=X2.5.,xmax=X97.5.), position=ggstance::position_dodgev(height=0.5)) +
        xlab("Proportion time Online > Telephone") +
        scale_color_manual(values=brewer.pal(3,"Set1")[-3],
                     guide=guide_legend(title=NULL)) +
  theme_bw() +
  theme(axis.text=element_text(color="black"),
        axis.title.y=element_blank(),
        axis.ticks=element_line(color="black"),
        panel.border=element_rect(color="black")) +
  geom_vline(xintercept=0.5, linetype="dashed")

png("output/state_scatter.png", width=7, height=7, units="in", res=300)
p
dev.off()

states = sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
prop_gt$ID = tolower(df$state[match(prop_gt$state,df$state.abb)])
states_merged = merge(states, prop_gt)

p = ggplot(states_merged, aes(fill=X50.)) +
        geom_sf(data = states, fill="lightgray") +    
        geom_sf() +
        scale_fill_gradient2(midpoint=0.5,
                             low = muted("purple"), mid = "white", high = muted("orange"),
                             guide=guide_colorbar(title="Proportion time Online > Telephone")) +
        facet_wrap(~ response, nrow=2) +
        theme_bw() +
        theme(axis.text=element_blank(),
              axis.ticks=element_blank(),
              legend.position="bottom")

png("output/state_map.png", width=5, height=6, units="in", res=300)
p
dev.off()

# I think the extremely wide credible intervals are because we are modeling each state, response, and mode
# combination separately.  Might be better to use a hierarchical model where we can borrow mode effect info.
# across states

prop_gt = data.frame(prop_gt, state_data[match(df$state[match(prop_gt$state,df$state.abb)],state_data$Geography),])
form = as.formula(paste0("X50. ~ ", paste(colnames(state_data[,-1]),collapse="+")))
mod = lm(form, data=prop_gt) # Online-Web decreases as HS.Grad increases (based on point estimates)
summary(mod)

summary(lm(X50. ~ median_age, data=prop_gt))

####################### 4. correlation analysis

# can make time series cross-correlation matrix for predictions associated with each mode...
# or for mode effect time series directly.  I'm using the latter approach for now.

# looks like interpolation/modeling is common used here but may not be ideal.
# Tools:
# https://rdrr.io/github/svdataman/sour/man/cross_correlate.html
# https://rdrr.io/cran/BINCOR/
# Papers/info:
# https://www.nonlin-processes-geophys.net/18/389/2011/npg-18-389-2011.pdf
# https://www.researchgate.net/publication/283076841_Correlation_analysis_for_unevenly_spaced_environmental_time_series
# https://stats.stackexchange.com/questions/7727/how-to-correlate-two-time-series-with-gaps-and-different-time-bases

# see also:
# https://stats.stackexchange.com/questions/133155/how-to-use-pearson-correlation-correctly-with-time-series
# maybe we can just diff and report Spearman's rho?

correlation_mat = expand.grid(response=as.character(unique(state_results$response)),
                              state1=as.character(unique(state_results$state)),
                              state2=as.character(unique(state_results$state)),
                              stringsAsFactors=FALSE) # ugh.. I hate factors

correlation_mat$rho = sapply(1:nrow(correlation_mat),function(i){
        print(i)                                     
        x = state_results$X50.[! is.na(state_results$time) &
                               state_results$response == correlation_mat$response[i] &
                               state_results$state == correlation_mat$state1[i]]
        y = state_results$X50.[! is.na(state_results$time) &
                               state_results$response == correlation_mat$response[i] &
                               state_results$state == correlation_mat$state2[i]]
        cor(diff(x),diff(y),method="spearman")
})

p = ggplot(correlation_mat, aes(x=state1,y=state2,fill=rho)) +
        geom_raster() +
        coord_equal() +
        facet_wrap(~ response, nrow=2) +
        theme_bw() +
        theme(axis.title=element_blank(),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
        scale_fill_gradient2(limits=c(-1,1),
                             low = muted("purple"), mid = "white", high = muted("orange"),
                             guide=guide_colorbar(title="Correlation"))

png("output/corr_mat_unsorted.png", width=5, height=7, units="in", res=300)
p
dev.off()

# we can perhaps shuffle the orders to make interpretation easier -- e.g.
# https://rdrr.io/cran/lessR/man/corReorder.html
# https://www.personality-project.org/r/html/mat.sort.html

plot_list = lapply(unique(correlation_mat$response), function(response){
        select_long = correlation_mat[correlation_mat$response == response,]
        select = dcast(select_long, state1 ~ state2, value.var="rho")
        rownames(select) = select$state1
        select$state1 = NULL
        ord = rownames(mat.sort(select)) # kind of weird to sort panels differently... maybe avoid
        #
        select_long$state1 = factor(select_long$state1, levels=ord)
        select_long$state2 = factor(select_long$state2, levels=ord)
        #
        p = ggplot(select_long, aes(x=state1,y=state2,fill=rho)) +
        geom_raster() +
        coord_equal() +
        theme_bw() +
        theme(axis.title=element_blank(),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
        scale_fill_gradient2(limits=c(-1,1),
                             low = muted("purple"), mid = "white", high = muted("orange"),
                             guide=guide_colorbar(title="Correlation"))
})

png("output/corr_mat.png", width=5, height=7, units="in", res=300)
grid.arrange(plot_list[[1]],plot_list[[2]],nrow=2)
dev.off()

####################### 5. prediction

# Maybe we can just use posterior distribution samples at the ends of the time series or something like that

