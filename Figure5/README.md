---
Title: The role of IAA and characterization of PIN transporters in complex streptophyte alga Chara braunii
Author: "Katarina Kurtović, Stanislav Vosolsobě, Daniel Nedvěd, Karel Müller, Petre Dobrev, Vojtěch Schmidt, Piotr Piszczek, Andre Kuhn, Adrijana Smoljan, Tom J. Fisher, Dolf Weijers, Jiří Friml, John L. Bowman & Jan Petrášek"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# ANALYSIS OF PIN FROM *CHARA BRAUNII* IN TOBACCO BY-2 CELLS
  
### R code by Stanislav Vosolsobě

## Loading a libraries

```r
library(peripheral)   # installed by http://kfrserver.natur.cuni.cz/lide/vosolsob/Peripheral/
library(emmeans)
```



## Fig. 5 b) Analysis of PM localisation of CbPINa and CbPINc

Detailed protocol was published in Vosolsobě, S., Schwarzerová, K., Petrášek, J. (2018): Determination of Plasma Membrane Partitioning for Peripherally-associated Proteins. JoVE (Journal of Visualized Experiments) (136), e57837

```r
profile <- read.profile(file="Results.csv")

# Now, the model of PM and cytoplasmic distribution based on previuos analysis is loaded
# For details, follow the protocol in http://kfrserver.natur.cuni.cz/lide/vosolsob/Peripheral/
load("model.RData")

filtered <- profile
md <- distrib(filtered, mp)
mdf <- filter(md, cut.off = 0.2)

df <- data.frame(rownames(mdf$affinity), mdf$affinity)
colnames(df) <- c(" ", colnames(mdf$affinity))
write.csv(df, file = "distrib.csv", row.names = FALSE, quote = FALSE)

pdf("distrib.pdf", width=10, height = 8)
boxplot(mdf, variants=c("var2"))
dev.off()

anova(lm(log(mdf$affinity$PM.aff)~mdf$affinity$var2))
# p = 3.898e-13 ***
```


