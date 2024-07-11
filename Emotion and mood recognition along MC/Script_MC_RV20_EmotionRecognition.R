pacman::p_load(psych, nlme, lme4, foreign, ggplot2, effects, TOSTER, multcomp, brms, mediation, Hmisc, BayesFactor, smooth, broom, lsmeans, forecast, zoo, dplyr, simr, longpower)

setwd("D:/Eigene Dateien/Forschung/Manuskripte/submitted/$RV20_MC_Mood_Pletzer/")

########################################################################################################
####################################### HORMONES #######################################################
########################################################################################################

data_horm <- read.spss(file = "RV20_Questionnaires_MC_20230519.sav", to.data.frame = T,  use.value.labels = F)
data_horm <- subset(data_horm, include == 1)
data_horm$Phase_factor <- factor(data_horm$Phase)

data_horm_long <- subset(data_horm, VPN > 400 & VPN < 500) 
data_horm_cross <- subset(data_horm, Session == 1) 

#cycle data
describe.by(data_horm_long$cycle_day_forward, data_horm_long$Phase_factor) 
describe.by(data_horm_long$cycle_day_backward, data_horm_long$Phase_factor) 

describe.by(data_horm_cross$cycle_day_forward, data_horm_cross$Phase_factor) 
describe.by(data_horm_cross$cycle_day_backward, data_horm_cross$Phase_factor) 

#estradiol
Estradiol_long <- lme(data = data_horm_long,  method="REML", scale(estradiol) ~ 1, random = ~ 1|VPN) 
Estradiol_long <- update(Estradiol_long, .~.+Phase_factor)
anova(Estradiol_long, type = "marginal")
summary(glht(Estradiol_long, alternative = "two.sided", linfct = mcp(Phase_factor = "Tukey")))

describe.by(data_horm_long$estradiol, data_horm_long$Phase_factor) 

ggplot(data = data_horm_long, aes(x=Phase_factor, y=estradiol)) + 
  geom_boxplot(fill = "tomato1") + theme_bw(base_size = 20) + scale_y_continuous(limits = c(0,3), name = "estradiol") + scale_x_discrete(labels = c("menses","pre-ovulatory","luteal")) + theme(axis.text=element_text(size=20))

Estradiol_cross <- lm(data = data_horm_cross,  scale(estradiol) ~ Phase_factor) 
anova(Estradiol_cross)
summary(glht(Estradiol_cross, alternative = "two.sided", linfct = mcp(Phase_factor = "Tukey")))

describe.by(data_horm_cross$estradiol, data_horm_cross$Phase_factor) 

#progesterone
Prog_long <- lme(data = data_horm_long,  method="REML", scale(progesterone) ~ 1, random = ~ 1|VPN) 
Prog_long <- update(Prog_long, .~.+Phase_factor)
anova(Prog_long, type = "marginal")
summary(glht(Prog_long, alternative = "two.sided", linfct = mcp(Phase_factor = "Tukey")))

describe.by(data_horm_long$progesterone, data_horm_long$Phase_factor) 

ggplot(data = data_horm_long, aes(x=Phase_factor, y=progesterone)) + 
  geom_boxplot(fill = "dodgerblue") + theme_bw(base_size = 20) + scale_y_continuous(limits = c(0,500), name = "progesterone") + scale_x_discrete(labels = c("menses","pre-ovulatory","luteal")) + theme(axis.text=element_text(size=20))

Prog_cross <- lm(data = data_horm_cross,  scale(progesterone) ~ Phase_factor) 
anova(Prog_cross)
summary(glht(Prog_cross, alternative = "two.sided", linfct = mcp(Phase_factor = "Tukey")))

describe.by(data_horm_cross$progesterone, data_horm_cross$Phase_factor) 

#####################################################################################################################
######################################### Emotion Recognition #######################################################
#####################################################################################################################

data_emotion <- read.spss(file = "RV20_EmotionRecognition_MC_20230519.sav", to.data.frame = T,  use.value.labels = F)
data_emotion <- subset(data_emotion, include == 1)
data_emotion$Phase_factor <- factor(data_emotion$Phase)
data_emotion$Session_factor <- factor(data_emotion$Session)
data_emotion$Emotion_factor <- factor(data_emotion$EmotionNr2)

data_emotion_long <- subset(data_emotion, VPN > 400 & VPN < 500)
data_emotion_cross <- subset(data_emotion, Session == 1)

#outlier correction
cutoff_long = mean(data_emotion_long$RT) + 3*sd(data_emotion_long$RT)
data_emotion_long_noOutlier <- subset(data_emotion_long, RT < cutoff_long)

cutoff_cross = mean(data_emotion_cross$RT) + 3*sd(data_emotion_cross$RT)
data_emotion_cross_noOutlier <- subset(data_emotion_cross, RT < cutoff_cross)

################################################################################
############################# Menstrual cycle effects ##########################
################################################################################

############################################ Speed ##################################################################

data_emotion_long_corr <- subset(data_emotion_long_noOutlier, Hit == 1)
MergedData_long_RT <- aggregate(RT ~ VPN + Session + Phase_factor + EmotionNr + Emotion_factor + progesterone + estradiol, FUN = "mean", data = data_emotion_long_corr)
MergedData_long_RT$VPN_factor <- factor(MergedData_long_RT$VPN)

data_emotion_cross_corr <- subset(data_emotion_cross_noOutlier, Hit == 1)
MergedData_cross_RT <- aggregate(RT ~ VPN + Session + Phase_factor + EmotionNr + Emotion_factor + progesterone + estradiol, FUN = "mean", data = data_emotion_cross_corr)
MergedData_cross_RT$VPN_factor <- factor(MergedData_cross_RT$VPN)

#cycle phase comparison
Emotion_RT_long <- lme(data = MergedData_long_RT,  method="REML", scale(RT) ~ Phase_factor*Emotion_factor+scale(Session), random = ~ 1|VPN) 
anova(Emotion_RT_long, type = "marginal") #no interaction session*cycle, no cycle effect
summary(glht(Emotion_RT_long, alternative = "two.sided", linfct = mcp(Emotion_factor = "Tukey")))

bfWithInteraction_long <- lmBF(RT ~ Phase_factor*Emotion_factor + Session + VPN_factor, data = MergedData_long_RT, whichRandom = "VPN") 
bfWithInteraction_long = recompute(bfWithInteraction_long, iterations = 10000)
bfWithPhase_long <- lmBF(RT ~ Phase_factor+Emotion_factor + Session + VPN_factor, data = MergedData_long_RT, whichRandom = "VPN") 
bfWithPhase_long = recompute(bfWithPhase_long, iterations = 10000)
bfWithoutPhase_long <- lmBF(RT ~ Emotion_factor + Session + VPN_factor, data = MergedData_long_RT, whichRandom = "VPN") 
bfWithoutPhase_long = recompute(bfWithoutPhase_long, iterations = 10000)

bfWithPhase_long/bfWithInteraction_long
bfWithoutPhase_long/bfWithPhase_long

ggplot(data = MergedData_long_RT, aes(x=Phase_factor, y=RT, fill = Emotion_factor)) + 
  geom_boxplot() + theme_bw(base_size = 20) + scale_y_continuous(name = "RT [ms]") + 
  scale_x_discrete(labels = c("menses","pre-ovulatory","luteal")) + theme(axis.text=element_text(size=20)) + 
  scale_fill_manual(values = c("white", "yellow", "grey", "red", "green", "blue"), labels = c("neutral", "happy", "fear", "angry", "disgust", "sad"))

Emotion_RT_cross <- lme(data = MergedData_cross_RT,  method="REML", scale(RT) ~ Phase_factor*Emotion_factor, random = ~ 1|VPN) 
anova(Emotion_RT_cross, type = "marginal") #no interaction session*cycle, no cycle effect
summary(glht(Emotion_RT_cross, alternative = "two.sided", linfct = mcp(Emotion_factor = "Tukey")))

bfWithInteraction_cross <- lmBF(RT ~ Phase_factor*Emotion_factor + VPN_factor, data = MergedData_cross_RT, whichRandom = "VPN") 
bfWithInteraction_cross = recompute(bfWithInteraction_cross, iterations = 10000)
bfWithPhase_cross <- lmBF(RT ~ Phase_factor+Emotion_factor + VPN_factor, data = MergedData_cross_RT, whichRandom = "VPN") 
bfWithPhase_cross = recompute(bfWithPhase_cross, iterations = 10000)
bfWithoutPhase_cross <- lmBF(RT ~ Emotion_factor + VPN_factor, data = MergedData_cross_RT, whichRandom = "VPN") 
bfWithoutPhase_cross = recompute(bfWithoutPhase_cross, iterations = 10000)

bfWithPhase_cross/bfWithInteraction_cross
bfWithoutPhase_cross/bfWithPhase_cross

#association to sex hormones
Emotion_RT_horm <- lme(data = MergedData_long_RT,  method="REML", scale(RT) ~ scale(estradiol)*scale(progesterone)*Emotion_factor+scale(Session), random = ~ 1|VPN) 
anova(Emotion_RT_horm, type = "marginal") 

######################################## Accuracy ##################################################################

MergedData_long_acc <- aggregate(Hit ~ VPN + Session + Phase_factor + EmotionNr + Emotion_factor + progesterone + estradiol, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_acc <- subset(MergedData_long_acc, Hit >= 2)
MergedData_long_acc$Accuracy <- MergedData_long_acc$Hit/10*100 
MergedData_long_acc$VPN_factor <- factor(MergedData_long_acc$VPN)

MergedData_cross_acc <- aggregate(Hit ~ VPN + Session + Phase_factor + EmotionNr + Emotion_factor + progesterone + estradiol, FUN = "sum", data = data_emotion_cross_noOutlier)
MergedData_cross_acc <- subset(MergedData_cross_acc, Hit >= 2)
MergedData_cross_acc$Accuracy <- MergedData_cross_acc$Hit/10*100 
MergedData_cross_acc$VPN_factor <- factor(MergedData_cross_acc$VPN)

#cycle phase comparison
Emotion_Acc_long <- lme(data = MergedData_long_acc,  method="REML", scale(Accuracy) ~ Phase_factor*Emotion_factor+scale(Session), random = ~ 1|VPN) 
anova(Emotion_Acc_long, type = "marginal") #no interaction session*cycle, no cycle effect
summary(glht(Emotion_Acc_long, alternative = "two.sided", linfct = mcp(Emotion_factor = "Tukey")))

bfWithInteraction_long <- lmBF(Accuracy ~ Phase_factor*Emotion_factor + Session + VPN_factor, data = MergedData_long_acc, whichRandom = "VPN") 
bfWithInteraction_long = recompute(bfWithInteraction_long, iterations = 10000)
bfWithPhase_long <- lmBF(Accuracy ~ Phase_factor+Emotion_factor + Session + VPN_factor, data = MergedData_long_acc, whichRandom = "VPN") 
bfWithPhase_long = recompute(bfWithPhase_long, iterations = 10000)
bfWithoutPhase_long <- lmBF(Accuracy ~ Emotion_factor + Session + VPN_factor, data = MergedData_long_acc, whichRandom = "VPN") 
bfWithoutPhase_long = recompute(bfWithoutPhase_long, iterations = 10000)

bfWithPhase_long/bfWithInteraction_long
bfWithoutPhase_long/bfWithPhase_long

ggplot(data = MergedData_long_acc, aes(x=Phase_factor, y=Accuracy, fill = Emotion_factor)) + 
  geom_boxplot() + theme_bw(base_size = 20) + scale_y_continuous(name = "Accuracy [%]") + 
  scale_x_discrete(labels = c("menses","pre-ovulatory","luteal")) + theme(axis.text=element_text(size=20)) + 
  scale_fill_manual(values = c("white", "yellow", "grey", "red", "green", "blue"), labels = c("neutral", "happy", "fear", "angry", "disgust", "sad"))

Emotion_Acc_cross <- lme(data = MergedData_cross_acc,  method="REML", scale(Accuracy) ~ Phase_factor*Emotion_factor, random = ~ 1|VPN) 
anova(Emotion_Acc_cross, type = "marginal") #no interaction session*cycle, no cycle effect
summary(glht(Emotion_Acc_cross, alternative = "two.sided", linfct = mcp(Emotion_factor = "Tukey")))

bfWithInteraction_cross <- lmBF(Accuracy ~ Phase_factor*Emotion_factor + VPN_factor, data = MergedData_cross_acc, whichRandom = "VPN") 
bfWithInteraction_cross = recompute(bfWithInteraction_cross, iterations = 10000)
bfWithPhase_cross <- lmBF(Accuracy ~ Phase_factor+Emotion_factor + VPN_factor, data = MergedData_cross_acc, whichRandom = "VPN") 
bfWithPhase_cross = recompute(bfWithPhase_cross, iterations = 10000)
bfWithoutPhase_cross <- lmBF(Accuracy ~ Emotion_factor + VPN_factor, data = MergedData_cross_acc, whichRandom = "VPN") 
bfWithoutPhase_cross = recompute(bfWithoutPhase_cross, iterations = 10000)

bfWithPhase_cross/bfWithInteraction_cross
bfWithoutPhase_cross/bfWithPhase_cross

#association to sex hormones
Emotion_Acc_horm_long <- lme(data = MergedData_long_acc,  method="REML", scale(Accuracy) ~ scale(estradiol)*scale(progesterone)*Emotion_factor+scale(Session), random = ~ 1|VPN) 
anova(Emotion_Acc_horm_long, type = "marginal") 

Emotion_Acc_horm_cross <- lme(data = MergedData_cross_acc,  method="REML", scale(Accuracy) ~ scale(estradiol)*scale(progesterone)*Emotion_factor, random = ~ 1|VPN) 
anova(Emotion_Acc_horm_cross, type = "marginal") 

############################################ Frequency ##############################################################

MergedData_long_freq_neutral <- aggregate(neutral ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_happy <- aggregate(happy ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_sad <- aggregate(sad ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_angry <- aggregate(angry ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_fear <- aggregate(fear ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_disgust <- aggregate(disgust ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol, FUN = "sum", data = data_emotion_long_noOutlier)

MergedData_long_freq_neutral <- rename(MergedData_long_freq_neutral, frequency = neutral)
MergedData_long_freq_happy <- rename(MergedData_long_freq_happy, frequency = happy)
MergedData_long_freq_sad <- rename(MergedData_long_freq_sad, frequency  = sad)
MergedData_long_freq_angry <- rename(MergedData_long_freq_angry, frequency = angry)
MergedData_long_freq_fear <- rename(MergedData_long_freq_fear, frequency = fear)
MergedData_long_freq_disgust <- rename(MergedData_long_freq_disgust, frequency = disgust)

MergedData_long_freq_neutral <- subset(MergedData_long_freq_neutral, frequency != 0)
MergedData_long_freq_happy <- subset(MergedData_long_freq_happy, frequency != 0)
MergedData_long_freq_sad <- subset(MergedData_long_freq_sad, frequency != 0)
MergedData_long_freq_angry <- subset(MergedData_long_freq_angry, frequency != 0)
MergedData_long_freq_fear <- subset(MergedData_long_freq_fear, frequency != 0)
MergedData_long_freq_disgust <- subset(MergedData_long_freq_disgust, frequency != 0)

MergedData_long_freq <- do.call("rbind", list(MergedData_long_freq_neutral, MergedData_long_freq_happy, MergedData_long_freq_sad, MergedData_long_freq_angry, MergedData_long_freq_fear, MergedData_long_freq_disgust))
MergedData_long_freq$Emotion_factor <- factor(MergedData_long_freq$Answer)
MergedData_long_freq$VPN_factor <- factor(MergedData_long_freq$VPN)

#cycle phase comparison
Emotion_freq_long <- lme(data = MergedData_long_freq,  method="REML", scale(frequency) ~ Phase_factor*Emotion_factor+scale(Session), random = ~ 1|VPN) 
anova(Emotion_freq_long, type = "marginal") #no interaction session*cycle, no cycle effect
summary(glht(Emotion_freq_long, alternative = "two.sided", linfct = mcp(Emotion_factor = "Tukey")))

bfWithInteraction <- lmBF(frequency ~ Phase_factor*Emotion_factor + Session + VPN_factor, data = MergedData_long_freq, whichRandom = "VPN") 
bfWithInteraction = recompute(bfWithInteraction, iterations = 10000)
bfWithPhase <- lmBF(frequency ~ Phase_factor+Emotion_factor + Session + VPN_factor, data = MergedData_long_freq, whichRandom = "VPN") 
bfWithPhase = recompute(bfWithPhase, iterations = 10000)
bfWithoutPhase <- lmBF(frequency ~ Emotion_factor + Session + VPN_factor, data = MergedData_long_freq, whichRandom = "VPN") 
bfWithoutPhase = recompute(bfWithoutPhase, iterations = 10000)

bfWithPhase/bfWithInteraction
bfWithoutPhase/bfWithPhase

ggplot(data = MergedData_long_freq, aes(x=Phase_factor, y=frequency, fill = Emotion_factor)) + 
  geom_boxplot() + theme_bw(base_size = 20) + scale_y_continuous(name = "frequency") + 
  scale_x_discrete(labels = c("menses","pre-ovulatory","luteal")) + theme(axis.text=element_text(size=20)) + 
  scale_fill_manual(values = c("red", "green", "grey", "yellow", "white", "blue"))

#association to sex hormones
Emotion_freq_horm <- lme(data = MergedData_long_freq,  method="REML", scale(frequency) ~ scale(estradiol)*scale(progesterone)*Emotion_factor+scale(Session), random = ~ 1|VPN) 
anova(Emotion_freq_horm, type = "marginal") 

################################################################################
############################# Exploratory: Moderation by mood ##################################
################################################################################

######################################## Speed ##################################################################

#moderation of overall emotion recognition by PMS
MergedData_long_RT_mood <- aggregate(RT ~ VPN + Session + Phase_factor + EmotionNr + Emotion_factor + progesterone + estradiol + PMS + Affect_Diff + STAI + DRSP4 + BAI_RW + BDI_RW, FUN = "mean", data = data_emotion_long_corr)
Emotion_RT_long_PMS <- lme(data = MergedData_long_RT_mood,  method="REML", scale(RT) ~ Phase_factor*Emotion_factor*scale(PMS)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_RT_long_PMS, type = "marginal") 

MergedData_cross_RT_mood <- aggregate(RT ~ VPN + Phase_factor + EmotionNr + Emotion_factor + progesterone + estradiol + PMS + Affect_Diff + STAI + DRSP4 + BAI_RW + BDI_RW, FUN = "mean", data = data_emotion_cross_corr)
Emotion_RT_cross_PMS <- lme(data = MergedData_cross_RT_mood,  method="REML", scale(RT) ~ Phase_factor*Emotion_factor*scale(PMS), random = ~ 1|VPN) 
anova(Emotion_RT_cross_PMS, type = "marginal") 

#moderation of sadness recognition by affect
MergedData_long_RT_sad <- subset(MergedData_long_RT_mood, EmotionNr == 2)
Emotion_RT_long_Affect_sad <- lme(data = MergedData_long_RT_sad,  method="REML", scale(RT) ~ Phase_factor*scale(Affect_Diff)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_RT_long_Affect_sad, type = "marginal") 

MergedData_cross_RT_sad <- subset(MergedData_cross_RT_mood, EmotionNr == 2)
Emotion_RT_cross_Affect_sad <- lme(data = MergedData_cross_RT_sad,  method="REML", scale(RT) ~ Phase_factor*scale(Affect_Diff), random = ~ 1|VPN) 
anova(Emotion_RT_cross_Affect_sad, type = "marginal") 

#moderation of fear recognition by state anxiety
MergedData_long_RT_fear <- subset(MergedData_long_RT_mood, EmotionNr == 3)
Emotion_RT_long_STAI_fear <- lme(data = MergedData_long_RT_fear,  method="REML", scale(RT) ~ Phase_factor*scale(STAI)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_RT_long_STAI_fear, type = "marginal") 

MergedData_cross_RT_fear <- subset(MergedData_cross_RT_mood, EmotionNr == 3)
Emotion_RT_cross_STAI_fear <- lme(data = MergedData_cross_RT_fear,  method="REML", scale(RT) ~ Phase_factor*scale(STAI), random = ~ 1|VPN) 
anova(Emotion_RT_cross_STAI_fear, type = "marginal") 

#moderation of anger recognition by irritability
MergedData_long_RT_angry <- subset(MergedData_long_RT_mood, EmotionNr == 4)
Emotion_RT_long_Irritability_angry <- lme(data = MergedData_long_RT_angry,  method="REML", scale(RT) ~ Phase_factor*scale(DRSP4)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_RT_long_Irritability_angry, type = "marginal") 

MergedData_cross_RT_angry <- subset(MergedData_cross_RT_mood, EmotionNr == 4)
Emotion_RT_cross_Irritability_angry <- lme(data = MergedData_cross_RT_angry,  method="REML", scale(RT) ~ Phase_factor*scale(DRSP4), random = ~ 1|VPN) 
anova(Emotion_RT_cross_Irritability_angry, type = "marginal") 

######################################## Accuracy ##################################################################

MergedData_long_acc_mood <- aggregate(Hit ~ VPN + Session + Phase_factor + EmotionNr + Emotion_factor + progesterone + estradiol + Affect_Diff + STAI + DRSP4 + PMS + BDI_RW + BAI_RW, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_acc_mood <- subset(MergedData_long_acc_mood, Hit >= 2 & Emotion_factor != 0 & Emotion_factor != 1)
MergedData_long_acc_mood$Accuracy <- MergedData_long_acc_mood$Hit/10*100 

MergedData_cross_acc_mood <- aggregate(Hit ~ VPN + Session + Phase_factor + EmotionNr + Emotion_factor + progesterone + estradiol + Affect_Diff + STAI + DRSP4 + PMS + BDI_RW + BAI_RW, FUN = "sum", data = data_emotion_cross_noOutlier)
MergedData_cross_acc_mood <- subset(MergedData_cross_acc_mood, Hit >= 2 & Emotion_factor != 0 & Emotion_factor != 1)
MergedData_cross_acc_mood$Accuracy <- MergedData_cross_acc_mood$Hit/10*100 

#moderation of overall emotion recognition by PMS
Emotion_Acc_long_PMS <- lme(data = MergedData_long_acc_mood,  method="REML", scale(Accuracy) ~ Phase_factor*Emotion_factor*scale(PMS)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_Acc_long_PMS, type = "marginal") #no interaction session*cycle, no cycle effect

Emotion_acc_cross_PMS <- lme(data = MergedData_cross_acc_mood,  method="REML", scale(Accuracy) ~ Phase_factor*Emotion_factor*scale(PMS), random = ~ 1|VPN) 
anova(Emotion_acc_cross_PMS, type = "marginal") 

#moderation of sadness recognition by affect
MergedData_long_acc_sad <- subset(MergedData_long_acc_mood, EmotionNr == 2)
Emotion_Acc_long_affect_sad <- lme(data = MergedData_long_acc_sad,  method="REML", scale(Accuracy) ~ Phase_factor*scale(Affect_Diff)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_Acc_long_affect_sad, type = "marginal") 

MergedData_cross_acc_sad <- subset(MergedData_cross_acc_mood, EmotionNr == 2)
Emotion_Acc_cross_Affect_sad <- lme(data = MergedData_cross_acc_sad,  method="REML", scale(Accuracy) ~ Phase_factor*scale(Affect_Diff), random = ~ 1|VPN) 
anova(Emotion_Acc_cross_Affect_sad, type = "marginal") 

#moderation of fear recognition by state anxiety
MergedData_long_acc_fear <- subset(MergedData_long_acc_mood, EmotionNr == 3)
Emotion_Acc_long_STAI_fear <- lme(data = MergedData_long_acc_fear,  method="REML", scale(Accuracy) ~ Phase_factor*scale(STAI)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_Acc_long_STAI_fear, type = "marginal") 

MergedData_cross_acc_fear <- subset(MergedData_cross_acc_mood, EmotionNr == 3)
Emotion_Acc_cross_STAI_fear <- lme(data = MergedData_cross_acc_fear,  method="REML", scale(Accuracy) ~ Phase_factor*scale(STAI), random = ~ 1|VPN) 
anova(Emotion_Acc_cross_STAI_fear, type = "marginal") 

#moderation of anger recognition by irritability
MergedData_long_acc_angry <- subset(MergedData_long_acc_mood, EmotionNr == 4)
Emotion_Acc_long_Irritability_angry <- lme(data = MergedData_long_acc_angry,  method="REML", scale(Accuracy) ~ Phase_factor*scale(DRSP4)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_Acc_long_Irritability_angry, type = "marginal") 

MergedData_cross_acc_angry <- subset(MergedData_cross_acc_mood, EmotionNr == 4)
Emotion_Acc_cross_Irritability_angry <- lme(data = MergedData_cross_acc_angry,  method="REML", scale(Accuracy) ~ Phase_factor*scale(DRSP4), random = ~ 1|VPN) 
anova(Emotion_Acc_cross_Irritability_angry, type = "marginal") 

################################# Frequency ##########################################################################

MergedData_long_freq_neutral_mood <- aggregate(neutral ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol + PMS + Affect_Diff + STAI + DRSP4 + BAI_RW + BDI_RW, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_happy_mood <- aggregate(happy ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol + PMS + Affect_Diff + STAI + DRSP4 + BAI_RW + BDI_RW, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_sad_mood <- aggregate(sad ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol + PMS + Affect_Diff + STAI + DRSP4 + BAI_RW + BDI_RW, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_angry_mood <- aggregate(angry ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol + PMS + Affect_Diff + STAI + DRSP4 + BAI_RW + BDI_RW, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_fear_mood <- aggregate(fear ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol + PMS + Affect_Diff + STAI + DRSP4 + BAI_RW + BDI_RW, FUN = "sum", data = data_emotion_long_noOutlier)
MergedData_long_freq_disgust_mood <- aggregate(disgust ~ VPN + Session + Phase_factor + Answer + progesterone + estradiol + PMS + Affect_Diff + STAI + DRSP4 + BAI_RW + BDI_RW, FUN = "sum", data = data_emotion_long_noOutlier)

MergedData_long_freq_neutral_mood <- rename(MergedData_long_freq_neutral_mood, frequency = neutral)
MergedData_long_freq_happy_mood <- rename(MergedData_long_freq_happy_mood, frequency = happy)
MergedData_long_freq_sad_mood <- rename(MergedData_long_freq_sad_mood, frequency  = sad)
MergedData_long_freq_angry_mood <- rename(MergedData_long_freq_angry_mood, frequency = angry)
MergedData_long_freq_fear_mood <- rename(MergedData_long_freq_fear_mood, frequency = fear)
MergedData_long_freq_disgust_mood <- rename(MergedData_long_freq_disgust_mood, frequency = disgust)

MergedData_long_freq_neutral_mood <- subset(MergedData_long_freq_neutral_mood, frequency != 0)
MergedData_long_freq_happy_mood <- subset(MergedData_long_freq_happy_mood, frequency != 0)
MergedData_long_freq_sad_mood <- subset(MergedData_long_freq_sad_mood, frequency != 0)
MergedData_long_freq_angry_mood <- subset(MergedData_long_freq_angry_mood, frequency != 0)
MergedData_long_freq_fear_mood <- subset(MergedData_long_freq_fear_mood, frequency != 0)
MergedData_long_freq_disgust_mood <- subset(MergedData_long_freq_disgust_mood, frequency != 0)

MergedData_long_freq_mood <- do.call("rbind", list(MergedData_long_freq_neutral_mood, MergedData_long_freq_happy_mood, MergedData_long_freq_sad_mood, MergedData_long_freq_angry_mood, MergedData_long_freq_fear_mood, MergedData_long_freq_disgust_mood))
MergedData_long_freq_mood$Emotion_factor <- factor(MergedData_long_freq_mood$Answer)
MergedData_long_freq_mood$VPN_factor <- factor(MergedData_long_freq_mood$VPN)

#moderation of overall emotion recognition by PMS
Emotion_freq_long_PMS <- lme(data = MergedData_long_freq_mood,  method="REML", scale(frequency) ~ Phase_factor*Emotion_factor*scale(PMS)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_freq_long_PMS, type = "marginal") 

#moderation of sadness recognition by affect
Emotion_freq_long_Affect_sad <- lme(data = MergedData_long_freq_sad_mood,  method="REML", scale(frequency) ~ Phase_factor*scale(Affect_Diff)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_freq_long_Affect_sad, type = "marginal") 

#moderation of fear recognition by state anxiety
Emotion_freq_long_STAI_fear <- lme(data = MergedData_long_freq_fear_mood,  method="REML", scale(frequency) ~ Phase_factor*scale(STAI)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_freq_long_STAI_fear, type = "marginal") 

#moderation of anger recognition by irritability
Emotion_freq_long_Irritability_angry <- lme(data = MergedData_long_freq_angry_mood,  method="REML", scale(frequency) ~ Phase_factor*scale(DRSP4)+scale(Session), random = ~ 1|VPN) 
anova(Emotion_freq_long_Irritability_angry, type = "marginal") 

#####################################################################################################################
######################################### MOOD ######################################################################
#####################################################################################################################

data_mood <- read.spss(file = "RV20_Mood_MC_20230519.sav", to.data.frame = T,  use.value.labels = F)
data_mood <- subset(data_mood, include == 1)
data_mood$Phase_factor <- factor(data_mood$Phase)
data_mood$Symptom_type <- factor(data_mood$Valence, labels = c("psychological", "physical"))
data_mood$Valence_factor <- factor(data_mood$Valence, labels = c("positive", "negative"))

data_mood_long <- subset(data_mood, VPN > 400 & VPN < 500)

############################ premenstrual symptoms #####################################################

data_DRSP_long <- data_mood_long[!is.na(data_mood_long$DRSP_avg),]

#cycle phase comparison
DRSP_long <- lme(data = data_DRSP_long,  method="REML", scale(DRSP_avg) ~ 1, random = ~ 1|VPN) 
DRSP_long <- update(DRSP_long, .~.+Phase_factor*Symptom_type)
anova(DRSP_long)
summary(glht(DRSP_long, alternative = "two.sided", linfct = mcp(Phase_factor = "Tukey")))

ggplot(data = data_DRSP_long, aes(x=Phase_factor, y=DRSP_avg, fill = Symptom_type)) + 
  geom_boxplot() + theme_bw(base_size = 20) + scale_y_continuous(name = "premenstrual symptoms", limits = c(1,6)) + scale_x_discrete(labels = c("menses","pre-ovulatory","luteal")) + theme(axis.text=element_text(size=20))

#association with sex hormones
DRSP_horm <- lme(data = data_DRSP_long,  method="REML", scale(DRSP_avg) ~ 1, random = ~ 1|VPN) 
DRSP_horm <- update(DRSP_horm, .~.+scale(estradiol)*scale(progesterone)*Valence_factor)
summary(DRSP_horm)

#association with self-reported PMS (trait)
data_DRSP_long_PSST <- data_DRSP_long[!is.na(data_DRSP_long$PMS),]

DRSP_longitudinal_PSST <- lme(data = data_DRSP_long_PSST,  method="REML", scale(DRSP_avg) ~ 1, random = ~ 1|VPN) 
DRSP_longitudinal_PSST <- update(DRSP_longitudinal_PSST, .~.+Phase_factor*Valence_factor*scale(PMS))
anova(DRSP_longitudinal_PSST)

######################################### affect ####################################################################

data_affect_long <- data_mood_long[!is.na(data_mood_long$Affect),]

#cycle phase comparison
Affect_long <- lme(data = data_affect_long,  method="REML", scale(Affect) ~ 1, random = ~ 1|VPN) 
Affect_long <- update(Affect_long, .~.+Phase_factor*Valence_factor)
anova(Affect_long)

ggplot(data = data_affect_long, aes(x=Phase_factor, y=Affect, fill = Valence_factor)) + 
  geom_boxplot() + theme_bw(base_size = 20) + scale_y_continuous(name = "affect", limits = c(1,5.5)) + scale_x_discrete(labels = c("menses","pre-ovulatory","luteal")) + theme(axis.text=element_text(size=20))

#resolve interaction
data_affect_long_pos <- subset(data_affect_long, Valence == 0)
data_affect_long_neg <- subset(data_affect_long, Valence == 1)

Affect_long_pos <- lme(data = data_affect_long_pos,  method="REML", scale(Affect) ~ 1, random = ~ 1|VPN) 
Affect_long_pos <- update(Affect_long_pos, .~.+Phase_factor)
anova(Affect_long_pos)
summary(glht(Affect_long_pos, alternative = "two.sided", linfct = mcp(Phase_factor = "Tukey")))

Affect_long_neg <- lme(data = data_affect_long_neg,  method="REML", scale(Affect) ~ 1, random = ~ 1|VPN) 
Affect_long_neg <- update(Affect_long_neg, .~.+Phase_factor)
anova(Affect_long_neg)
summary(glht(Affect_long_neg, alternative = "two.sided", linfct = mcp(Phase_factor = "Tukey")))

#associations to sex hormones
Affect_horm_pos <- lme(data = data_affect_long_pos,  method="REML", scale(Affect) ~ 1, random = ~ 1|VPN) 
Affect_horm_pos <- update(Affect_horm_pos, .~.+scale(estradiol)*scale(progesterone))
summary(Affect_horm_pos) 

Affect_horm_neg <- lme(data = data_affect_long_neg,  method="REML", scale(Affect) ~ 1, random = ~ 1|VPN) 
Affect_horm_neg <- update(Affect_horm_neg, .~.+scale(estradiol)*scale(progesterone))
summary(Affect_horm_neg) 

#association to BDI
data_affect_long_BDI <- data_affect_long[!is.na(data_affect_long$BDI_RW),]
data_affect_long_BDI_pos <- subset(data_affect_long_BDI, Valence == 0)
data_affect_long_BDI_neg <- subset(data_affect_long_BDI, Valence == 1)

Affect_long_BDI_pos <- lme(data = data_affect_long_BDI_pos,  method="REML", scale(Affect) ~ 1, random = ~ 1|VPN) 
Affect_long_BDI_pos <- update(Affect_long_BDI_pos, .~.+Phase_factor*scale(BDI_RW))
anova(Affect_long_BDI_pos)
summary(Affect_long_BDI_pos)

Affect_long_BDI_neg <- lme(data = data_affect_long_BDI_neg,  method="REML", scale(Affect) ~ 1, random = ~ 1|VPN) 
Affect_long_BDI_neg <- update(Affect_long_BDI_neg, .~.+Phase_factor*scale(BDI_RW))
anova(Affect_long_BDI_neg)
summary(Affect_long_BDI_neg)

############################################## anxiety #############################################################

data_STAI_long <- data_mood_long[!is.na(data_mood_long$STAI),]

#cycle phase comparison
STAI_long <- lme(data = data_STAI_long,  method="REML", scale(STAI) ~ 1, random = ~ 1|VPN) 
STAI_long <- update(STAI_long, .~.+Phase_factor)
anova(STAI_long)
summary(glht(STAI_long, alternative = "two.sided", linfct = mcp(Phase_factor = "Tukey")))

ggplot(data = data_STAI_long, aes(x=Phase_factor, y=STAI)) + 
  geom_boxplot(fill = "dodgerblue") + theme_bw(base_size = 20) + scale_y_continuous(name = "anxiety") + scale_x_discrete(labels = c("menses","pre-ovulatory","luteal")) + theme(axis.text=element_text(size=20))

#association with sex hormones
STAI_horm <- lme(data = data_STAI_long,  method="REML", scale(STAI) ~ 1, random = ~ 1|VPN) 
STAI_horm <- update(STAI_horm, .~.+scale(estradiol)*scale(progesterone))
summary(STAI_horm) 

#modulation by trait anxiety
data_STAI_long_BAI <- data_STAI_long[!is.na(data_STAI_long$BAI_RW),]

STAI_longitudinal_BAI <- lme(data = data_STAI_long_BAI,  method="REML", scale(STAI) ~ 1, random = ~ 1|VPN) 
STAI_longitudinal_BAI <- update(STAI_longitudinal_BAI, .~.+Phase_factor*scale(BAI_RW))
anova(STAI_longitudinal_BAI)
summary(STAI_longitudinal_BAI)
