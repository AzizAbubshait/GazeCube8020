library(tidyverse)
library(lme4)
library(lmerTest)
library(afex)

# colors for the graph ----
cbbPalette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To do  ----
# read in data ----
myPath = paste(getwd(), "/Behavioral data iCube/", sep = "")
filenames = list.files(myPath)

all_dat = lapply(paste(myPath, filenames, sep = ""), read.csv)
all_dat = do.call(plyr::rbind.fill, all_dat)

all_dat %>%
  group_by(grouping) %>%
  summarize(n= n_distinct(participant))

data.frame(all_dat %>%
             group_by(participant, grouping) %>%
             summarize(n = n_distinct()))

# fix participant IDS
all_dat$participant[which(all_dat$participant==1092)] = 109
all_dat$participant[which(all_dat$participant==1282)] = 128
all_dat$participant[which(all_dat$participant==1491)] = 149
all_dat$participant[which(all_dat$participant==0)] = 110


p_dat = all_dat %>%
  filter(is.na(N) == F) %>%
  group_by(participant, gazeCon, grouping, diff) %>%
  #filter(dropResp.rt < 2) %>%
  mutate(error_rate = sum(dropResp.corr==0)/length(dropResp.corr)) %>%
  ungroup(gazeCon, grouping) %>%
  mutate(n_trial = row_number()) %>% ungroup(participant)

data.frame(p_dat %>%
             group_by(gazeCon, grouping, participant) %>%
             summarize(n_trial = n()))

p_dat %>% group_by(gazeCon, grouping) %>% summarize(n=n_distinct(participant))
length(unique(p_dat$participant))
unique(p_dat$participant)

# error rates ----
err_dat = p_dat %>% 
  mutate(diff_vis = ifelse(diff < 4, "High", ifelse(diff > 4, "Low", "NA"))) %>%
  group_by(participant, gazeCon, grouping, diff_vis) %>%
  summarize(error_rate = unique(error_rate)) %>%
    mutate(high = ifelse(error_rate > .3, "Hi", "Low"))
data.frame(err_dat)

err_dat %>%
  ggplot(aes(y=error_rate, grouping, color = gazeCon))+
  scale_color_manual(values = cbbPalette)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1, 
               position = position_dodge(.5))+
  stat_summary(fun.data = mean_se, geom = "point", size = 5,
               position = position_dodge(.5))+theme_bw()+
  facet_wrap(~diff_vis)+ggtitle("Plot 1")

lm_err = lmer(error_rate ~ grouping*gazeCon*diff_vis + (1|participant), err_dat)
anova(lm_err)
aov_err = aov_car(error_rate ~ grouping*gazeCon*diff_vis + Error(participant/gazeCon*diff_vis), err_dat)
summary(aov_err)

# RTs ----
p_dat %>%
  group_by(participant, grouping, gazeCon) %>%
  filter(dropResp.corr == 1,
         dropResp.rt < 10) %>%
  ggplot(aes(dropResp.rt)) +
  geom_histogram()

median(processed_dat2$dropResp.rt)
mean(processed_dat2$dropResp.rt)

processed_dat2 = p_dat %>%
  group_by(participant, grouping, gazeCon) %>%
  filter(dropResp.corr == 1,
         dropResp.rt < 2) %>%
  mutate(z_drop_within = scale(dropResp.rt)) %>%
  filter(abs(z_drop_within) < 4)

data.frame(processed_dat2 %>%
             group_by(gazeCon, grouping, participant) %>%
             summarize(n_trial = n()))

processed_dat2 %>%
  group_by(gazeCon, grouping, participant) %>%
  summarize(n_trial = n())

# plot trials without averaging:
ggplot(processed_dat2, aes(y = dropResp.rt, x = grouping, color = gazeCon, fill = gazeCon))+
  scale_fill_manual(values = cbbPalette)+
  scale_color_manual(values = cbbPalette)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1, 
               position = position_dodge(.5), color = "black")+
  stat_summary(fun.data = mean_se, geom = "point", size = 5,
               position = position_dodge(.5))+theme_bw()

lm1 = lmer(dropResp.rt ~ grouping*gazeCon + 
             (1|participant)+(1|n_trial), processed_dat2)
summary(lm1)

avg_dat = processed_dat2 %>%
  group_by(participant, gazeCon, grouping) %>%
  summarize(dropResp.rt = mean(dropResp.rt))

aov1 = aov_car(dropResp.rt ~ grouping*gazeCon + Error(participant/gazeCon),
               avg_dat)
summary(aov1)

# incorporate difficulty ----
processed_dat2 %>%
  mutate(diff_vis = ifelse(diff < 4, "High", ifelse(diff > 4, "Low", "NA"))) %>%
  ggplot(aes(y = dropResp.rt, x = grouping, color = gazeCon, fill = gazeCon))+
  scale_fill_manual(values = cbbPalette)+
  scale_color_manual(values = cbbPalette)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1, 
               position = position_dodge(.5), color = "black")+
  stat_summary(fun.data = mean_se, geom = "point", size = 5,
               position = position_dodge(.5))+theme_bw()+
  facet_wrap(~diff_vis)

lm2 = lmer(dropResp.rt ~ grouping*gazeCon*diff + 
             (1|participant)+(1|n_trial), processed_dat2)
anova(lm2)

processed_dat2 %>%
  mutate(diff_vis = ifelse(diff < 4, "High", ifelse(diff > 4, "Low", "NA"))) %>%
  ggplot(aes(y = dropResp.rt, x = diff_vis, color = gazeCon))+
  scale_color_manual(values = cbbPalette)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1, 
               position = position_dodge(.5))+
  stat_summary(fun.data = mean_se, geom = "point", size = 5,
               position = position_dodge(.5))+theme_bw()+
  ggtitle("Plot 2")

processed_dat2 %>%
  mutate(diff_vis = ifelse(diff < 4, "High", ifelse(diff > 4, "Low", "NA"))) %>%
  ggplot(aes(y = dropResp.rt, x = diff_vis, color = grouping))+
  scale_color_manual(values = cbbPalette)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1, 
               position = position_dodge(.5))+
  stat_summary(fun.data = mean_se, geom = "point", size = 5,
               position = position_dodge(.5))+theme_bw()+
  ggtitle("Plot 3")

processed_dat2 %>%
  mutate(diff_vis = ifelse(diff < 4, "High", ifelse(diff > 4, "Low", "NA"))) %>%
  ggplot(aes(y = dropResp.rt, x = grouping, color = gazeCon))+
  scale_color_manual(values = cbbPalette)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1, 
               position = position_dodge(.5))+
  stat_summary(fun.data = mean_se, geom = "point", size = 5,
               position = position_dodge(.5))+theme_bw()+
  facet_wrap(~diff_vis)+
  ggtitle("Overall Plot")

avg_dat_diff = processed_dat2 %>%
  mutate(diff_vis = ifelse(diff < 4, "High", ifelse(diff > 4, "Low", "NA"))) %>%
  group_by(participant, gazeCon, grouping,diff_vis) %>%
  summarize(dropResp.rt = mean(dropResp.rt))

aov2 = aov_car(dropResp.rt ~ grouping*gazeCon*diff_vis + 
                 Error(participant/gazeCon*diff_vis),
               avg_dat_diff)
summary(aov2)

