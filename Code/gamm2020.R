require(sf)
library(tidyverse)
library(scales)
library(forcats)
require(lubridate)
require(reshape2)
require(reshape)
require(boot)
require(mgcv)
require(ggplot2)
require(RColorBrewer)
require(plyr)
require(Hmisc)
require(gridExtra)
require(multcomp)
require(broom)
require(mgcv.helper) # devtools::install_github("samclifford/mgcv.helper")
require(tidycensus)
require(urbnmapr)
require(itsadug)

setwd("/home/chrisgraywolf/shared/analysis/misc/Heather/2021_conference/")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

####################### 0. clean up the election data
# https://stackoverflow.com/questions/46123285/printing-intermediate-results-without-breaking-pipeline-in-tidyverse

df_survey = read.csv("2020 Presidential Election Polls.csv", as.is=TRUE) %>%
        mutate(startdate = as.Date(format(as.Date(start_date, "%m/%d/%Y"), "20%y/%m/%d")),
               enddate = as.Date(format(as.Date(end_date, "%m/%d/%Y"), "20%y/%m/%d"))) %>%
        mutate(time=as.numeric(enddate)) %>%
        filter(methodology != "") %>%
        mutate(Mode = factor(str_detect(methodology,"Live"), # TODO - ask about Live Phone/Online/Text etc.
                             levels = c(FALSE, TRUE),
                             labels = c("Other","Live")),
               pollster = factor(pollster, levels=c(unique(pollster),"new")), # add "new" level for prediction
               state = gsub(" CD-.*","",state),
               state_abb = state.abb[match(state, state.name)]) %>%
        mutate(state_abb = case_when(state %in% "District of Columbia" ~ "DC",
                                     state %in% "" ~ "National",
                                     TRUE ~ state_abb)) %>%
        filter(answer %in% c("Biden","Trump")) %>%
        spread("answer", "pct") %>%
        dplyr::rename(adjpoll_trump = Trump,
                      adjpoll_biden = Biden)

# add some state-level covariates https://walker-data.com/tidycensus/articles/basic-usage.html
# https://www.socialexplorer.com/data/ACS2016_5yr/metadata/?ds=ACS16_5yr&var=B02001002
census_api_key("c2b4f05d5b5db357358c6d0105be651a02e8cc02")
acs_vars = load_variables(2019, "acs5", cache = TRUE)
# View(acs_vars)

census_vars = c("n_white"="B02001_002",
                "n_total"="B02001_001",
                "median_age"="B01002_001",
                "median_income"="B19326_001",
                "n_bach"="B15012_001",
                "n_male"="B01001_002",
                "n_foreign"="B05012_003",
                "n_total_pov"="B17001_001",
                "n_pov"="B17001_002") # TODO - consider other predictors

df_list = lapply(names(census_vars), function(var_name){
        print(var_name)
        out = get_acs(geography = "state", variables = census_vars[var_name], geometry = FALSE)
        out[,var_name] = out$estimate
        out = data.frame(out)
        return(out[,c("NAME",var_name)])
})

covar_df = reduce(df_list, inner_join, by = "NAME") %>%
        mutate(pct_white = 100 * n_white / n_total,
               pct_foreignborn = 100 * n_foreign / n_total, # TODO need to double check denominators
               pct_below_poverty = 100 * n_pov / n_total_pov,
               pct_bachorhigher = 100 * n_bach / n_total,
               sex_ratio_mw = n_male / (n_total - n_male),
               state_abb = state.abb[match(NAME, state.name)]) %>%
        filter(NAME != "Puerto Rico")
covar_df$state_abb[covar_df$NAME %in% "District of Columbia"] = "DC"

df = merge(df_survey, covar_df, by="state_abb", all.x=TRUE, all.y=FALSE)

# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/42MVDX
results = read.csv("1976-2020-president.csv", as.is=TRUE)
results = results[results$year == 2020 & results$candidate %in% c("BIDEN, JOSEPH R. JR","TRUMP, DONALD J.")
                  & ! results$writein & results$party_detailed %in% c("DEMOCRAT","REPUBLICAN"),]
results = spread(results[,c("state","candidate","candidatevotes")],"candidate","candidatevotes")
trump_states = results$state[results$`BIDEN, JOSEPH R. JR` < results$`TRUMP, DONALD J.`] # TODO -- account for split states

####################### 1. additive (across time) mode effect models (easier to estimate parameters)
# In short, for "Live" mode, all covars. but S(time, by=state) are zero, so we just treat
# this as the state-varying baseline.  Then, on top of that, we can add some "mode effect"
# for the "Other" mode, which can depend on covariates (e.g., age, income).
# ref.: https://rdrr.io/cran/mgcv/man/random.effects.html

df[df$Mode == "Live", c("median_age","median_income","pct_below_poverty","pct_bachorhigher",
                          "sex_ratio_mw","pct_foreignborn","pct_white")] = 0 # treat "Live" as baseline
df$is_web = df$Mode == "Other"
df$state = factor(df$state)

mod = bam(adjpoll_trump/1e2 ~ s(time,by=state) + is_web + median_age + median_income + pct_below_poverty +
          pct_bachorhigher + sex_ratio_mw + pct_foreignborn + pct_white + s(pollster,bs="re"),
          data=df[df$state != "",], # exclude national polls for this part
          family=betar(link="logit"),
          drop.unused.levels=FALSE,
          discrete=TRUE)

summary(mod)

par = c("is_webTRUE","median_age","median_income", "pct_below_poverty","pct_bachorhigher","sex_ratio_mw","pct_foreignborn","pct_white")
par_df =  mgcv.helper:::confint.gam(mod, parm=par) %>% data.frame
par_df[,-1] = round(par_df[,-1], 3) # TODO - consider back-transforming for interpretability
write.csv(par_df, "output/state_model_covar_results.csv", row.names=FALSE)

### now do state level mode effect predictions/CIs.
# so, for each row in the glht matrix, we want is_web to be TRUE and
# state-level covars. for a different state

# https://stat.ethz.ch/pipermail/r-help/2012-August/320895.html
test_mat = 0 * model.matrix(mod)[1:51,]
test_mat[,"is_webTRUE"] = 1 # intercept is always included in mode effect part
covars = intersect(colnames(test_mat), colnames(covar_df))
test_mat[,covars] = as.matrix(covar_df[,covars])

# simultaneous CIs: https://www.rdocumentation.org/packages/multcomp/versions/1.4-10/topics/glht-methods
conf_df = tidy(confint(glht(mod,linfct=test_mat))) %>% data.frame
conf_df$state = covar_df$NAME
conf_df$state = factor(conf_df$state, levels=conf_df$state[order(conf_df$estimate)])
conf_df$trump = toupper(conf_df$state) %in% trump_states # num. Trump votes > num. Biden votes ?
write.csv(conf_df, "output/conf_df.csv", row.names=FALSE)

p = ggplot(conf_df, aes(x=estimate,y=state,color=trump)) + # logit scale
        geom_point() +
        geom_errorbarh(aes(xmin=conf.low, xmax=conf.high)) +
        theme_bw() +
  xlab("Mode effect (Non-live - Live)") + # was: (Online - Telephone)
  scale_color_manual(values=brewer.pal(3,"Set1")[2:1], guide=FALSE) +
  theme(axis.text=element_text(color="black"),
        axis.title.y=element_blank(),
        axis.ticks=element_line(color="black"),
        panel.border=element_rect(color="black")) +
  geom_vline(xintercept=0,linetype="dashed") # shy Trump voters mostly in R states.. kind of weird

png("output/additive_mode.png", width=7, height=7, units="in", res=300)
p
dev.off()

pdf("output/additive_mode.pdf", width=7, height=7)
p
dev.off()

# make table of state level mode effect results (positive / negative / non-sig.)
conf_df$not_sig = conf_df$conf.low < 0 & conf_df$conf.high > 0
conf_df$outcome = NA
conf_df$outcome[conf_df$estimate < 0] = "Negative"
conf_df$outcome[conf_df$estimate > 0] = "Positive"
conf_df$outcome[conf_df$not_sig] = "Not Significant"

tab = table(ifelse(conf_df$trump,"Trump","Biden"), conf_df$outcome)
write.csv(tab, "output/state_results_tab.csv")

states <- get_urbn_map("states", sf = TRUE)
states_merged = merge(states, conf_df, by.x="state_name", by.y="state")

p = ggplot(states_merged[! states_merged$not_sig,], aes(fill=estimate)) +
        geom_sf(data=states, fill="lightgray") +
        geom_sf() +
        scale_fill_gradient2(midpoint=0,
                             low = muted("purple"), mid = "white", high = muted("orange"),
                             guide=guide_colorbar(title="Mode effect (Non-live - Live)")) +
        theme_bw() +
        theme(axis.text=element_blank(),
              axis.ticks=element_blank(),
              legend.position="bottom")

png("output/state_map.png", width=6, height=4, units="in", res=300)
p
dev.off()

pdf("output/state_map.pdf", width=6, height=4)
p
dev.off()

####################### 2. test for national-level mode effect w/ compareML

for(response in c("adjpoll_trump","adjpoll_biden")){ # national models
        print(response)
        df$y = df[,response]/1e2

        national = gam(y ~ s(time) + s(pollster,bs="re"),
                  family=betar(link="logit"),
                       data=df[df$state == "",], method="ML")
        national_mode = gam(y ~ s(time,by=Mode) + s(pollster,bs="re"), 
                family=betar(link="logit"),
                 data=df[df$state == "",], method="ML")
        compareML(national, national_mode)
}

####################### 3. predict and plot (national, with and without mode)

pred = do.call("rbind", lapply(c("adjpoll_trump","adjpoll_biden"), function(response){
        print(response)
        df$y = df[,response]/1e2

        national = gam(y ~ s(time) + s(pollster,bs="re"), data=df[df$state == "",],
                          family=betar(link="logit"),
                       method="ML",drop.unused.levels=FALSE)
        pred = data.frame(time = seq(min(df$time),max(df$time),length=1e2), Mode="Overall", pollster="new")
        pred$y = predict(national,newdata=pred,type="response")

        national_mode = gam(y ~ s(time,by=Mode) + s(pollster,bs="re"), data=df[df$state == "",],
                                  family=betar(link="logit"),
                            method="ML",drop.unused.levels=FALSE)
        pred_mode = data.frame(expand.grid(time = seq(min(df$time),max(df$time),length=1e2), Mode=unique(df$Mode)),
                                         pollster="new")
        pred_mode$y = predict(national_mode,newdata=pred_mode,type="response")
        
        out = rbind(pred,pred_mode)
        out$variable = response
        return(out)
}))

poll_df = melt(df, measure.vars=c("adjpoll_trump","adjpoll_biden"))
poll_df$Mode = factor(poll_df$Mode, levels=c("Overall","Other","Live"),
                              labels=c("Overall","Non-live","Live"))

pred$Mode = factor(pred$Mode, levels=c("Overall","Other","Live"),
                              labels=c("Overall","Non-live","Live"))
pred$variable = factor(pred$variable, levels=c("adjpoll_biden","adjpoll_trump"), labels=c("Biden","Trump"))
poll_df$variable = factor(poll_df$variable, levels=c("adjpoll_biden","adjpoll_trump"), labels=c("Biden","Trump"))

p = ggplot(pred, aes(x=as.Date(time,origin="1970-01-01"), y=value, color=Mode, size=Mode)) + 
        facet_wrap(~ variable) +
        theme_bw() +
        geom_point(data=poll_df[poll_df$state == "",], alpha=0.25,size=1,show.legend=FALSE) +
        geom_line(data=pred, aes(y=y*1e2,size=Mode)) +
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

####################### 4. compare predictions with results
# note - just working with state-level point estimates for now,
# but could consider state-level CIs, or national model predictions, etc.

mod_state = bam(adjpoll_trump/1e2 ~ s(time,by=state) + is_web*state + s(pollster,bs="re"),
          data=df[df$state != "",], # exclude national polls for this part
          family=betar(link="logit"),
          drop.unused.levels=FALSE,
          discrete=TRUE)

pred = expand.grid(state=covar_df$NAME, time=max(df$time), is_web=c(FALSE, TRUE), pollster="new")

pred$y = 100 * predict(mod_state,newdata=pred,type="response") # TODO - check extreme preds w/ "mod" (maybe need logit link or w/e)

pred = spread(pred, "is_web", "y") %>%
        dplyr::rename("Live" = "FALSE",
                      "Not live" = "TRUE")        

results$prop_trump = results$`TRUMP, DONALD J.` / (results$`TRUMP, DONALD J.` + results$`BIDEN, JOSEPH R. JR`)

pred$prop_trump = 100 * results$prop_trump[match(toupper(pred$state), results$state)]

pred[,c("Live","Not live")] = pred$prop_trump - pred[,c("Live","Not live")]

p = ggplot(pred, aes(x=Live, y=`Not live`)) +
        geom_point() +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        coord_fixed() +
        geom_abline(intercept=0,slope=1) +
        theme_bw() +
        ggtitle("Model residuals (state level) for prop. Trump") # TODO - figure out why nearly all are positive

png("output/residuals_plot.png", width=5.5, height=4.5, units="in", res=300)
p
dev.off()

pdf("output/residuals_plot.pdf", width=5.5, height=4.5)
p
dev.off()

####################### conceptual figure

set.seed(1)
spl_df = data.frame(predict(smooth.spline(x=seq(0,1,length=6),
                                          y=c(0    ,-0.25     ,     -1    ,-2,    -2  , -1.5),
                                          spar=0.2),
                            x=seq(0,1,length=1000)))

ran_ef = data.frame(do.call("cbind", lapply(rnorm(10), function(x) spl_df$y + 0.2*x)))
ran_ef[,5:10] = ran_ef[,5:10] + 2
ran_ef$x = spl_df$x
ran_ef = melt(ran_ef, id.vars="x")

spl_df$y2 = spl_df$y + 2
mode_loc = spl_df[700,]
spl_df = melt(spl_df, measure.vars= c("y","y2"))

p = ggplot(spl_df, aes(y=value, x=x, group=variable)) +
        geom_line() +
        geom_line(data=ran_ef, color="red", alpha=0.2) +
        theme_classic() +
        theme(axis.text=element_text(color="black"),
                axis.ticks=element_line(color="black")) +
        xlab("Date") +
        ylab("Proportion Trump (logit scale)") +
        #geom_vline(xintercept=1, color="darkgreen", linetype="dashed")
        scale_y_continuous(breaks=NULL) +
        scale_x_continuous(breaks=1, labels="Election", limits=c(0,1.1)) +
        geom_segment(color="blue", size=1, x=mode_loc$x, xend=mode_loc$x, y=mode_loc$y, yend=mode_loc$y2,
                  arrow = arrow(length = unit(0.5, "cm"))) +
        annotate("text", x=mode_loc$x, y=mode_loc$y+1, label="beta[0] + beta[1] * x[1][i] + ... + beta[p] * x[p][i]", parse=TRUE, hjust=-.05,
                color="blue", size=5) +
        geom_text(data=spl_df[1,], x=0.1, y=-0.5, label="Live mode") +
        geom_text(data=spl_df[1,], x=0.1, y=1.5, label="Non-live mode")

png("output/conceptual_logit.png", width=5.5, height=4.5, units="in", res=300)
p
dev.off()

pdf("output/conceptual_logit.pdf", width=5.5, height=4.5)
p
dev.off()

p = ggplot(spl_df, aes(y=exp(value), x=x, group=variable)) +
        geom_line() +
        geom_line(data=ran_ef, color="red", alpha=0.2) +
        theme_classic() +
        theme(axis.text=element_text(color="black"),
                axis.ticks=element_line(color="black")) +
        xlab("Date") +
        ylab("Odds of voting for Trump (logit scale)") +
        #geom_vline(xintercept=1, color="darkgreen", linetype="dashed")
        scale_y_continuous(breaks=NULL,) +
        scale_x_continuous(breaks=1, labels="Election", limits=c(0,1.1)) +
        geom_segment(color="blue", size=1, x=mode_loc$x, xend=mode_loc$x, y=exp(mode_loc$y), yend=exp(mode_loc$y2),
                  arrow = arrow(length = unit(0.5, "cm"))) +
        annotate("text", x=mode_loc$x, y=exp(mode_loc$y+1)+0.1, label="'×'*e^{beta[0] + beta[1] * x[1][i] + ... + beta[p] * x[p][i]}", parse=TRUE, hjust=-.05,
                color="blue", size=5) +
        geom_text(data=spl_df[1,], x=0.1, y=exp(-0.5), label="Live mode") +
        geom_text(data=spl_df[1,], x=0.1, y=exp(1.5), label="Non-live mode")

png("output/conceptual_odds.png", width=5.5, height=4.5, units="in", res=300)
p
dev.off()

pdf("output/conceptual_odds.pdf", width=5.5, height=4.5)
p
dev.off()

p = ggplot(spl_df, aes(y=inv.logit(value), x=x, group=variable)) +
        geom_line() +
        geom_line(data=ran_ef, color="red", alpha=0.2) +
        theme_classic() +
        theme(axis.text=element_text(color="black"),
                axis.ticks=element_line(color="black")) +
        xlab("Date") +
        ylab("Proportion Trump") +
        scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(breaks=1, labels="Election", limits=c(0,1.1)) +
        geom_segment(color="blue", size=1, x=mode_loc$x, xend=mode_loc$x, y=mode_loc$y, yend=mode_loc$y2,
                  arrow = arrow(length = unit(0.5, "cm"))) +
        geom_text(data=spl_df[1,], x=0.1, y=inv.logit(-0.5), label="Live mode") +
        geom_text(data=spl_df[1,], x=0.1, y=inv.logit(1.5), label="Non-live mode")

png("output/conceptual_prob.png", width=5.5, height=4.5, units="in", res=300)
p
dev.off()

pdf("output/conceptual_prob.pdf", width=5.5, height=4.5)
p
dev.off()
