library(tidyverse)
file <- "D:/OneDrive/Documents/Manuscripts/in preparation - MFM BRCA2 tracking/BRCA2 quantification/ExpMP2009_007 BRCA2 quantification.txt"
intensity_values <- read_delim(file)%>%
  mutate(Concentration=Concentration/300*stock_concentration)


stock_concentration	<- 418 #nM



summarized_intensity_values <- intensity_values %>%
  group_by(Concentration)%>%
  summarize(Intensity=mean(Mean))


p <- intensity_values %>%
ggplot(aes(x=Concentration,y=Mean))+geom_point()+  geom_smooth(method=lm)
p

m <- lm(Intensity ~Concentration , summarized_intensity_values);
m



