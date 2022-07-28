#load packages
library(tidyverse)
library(ggbeeswarm)

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#ffffff"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.margin = unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#ffffff",fill="#ffffff"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#c00000","#1F497D","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#6599D9","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}



#load data
folder <- "D:/OneDrive/Data2/220623 ExpMP2206_005 BRCA2 in situ quantification"
intensity_values <- read_delim(file.path(folder,"220623 calibration data.txt"))
nucleus_data <- read_csv(file.path(folder,"3Dsegmentation_output_mask.csv"))
foci_data <- read_csv(file.path(folder,"3Dsegmentation_output_Watershed_Foci.csv"))

#parameters
Nav <- 6.02214076E23
Voxel_volume <- 0.1317882*0.1317882*0.4564125 #um3
stock_concentration	<- 418	#nM

#calculate concentrations of the Halotag complex in media 
intensity_values <- intensity_values %>%
  mutate(Concentration=Concentration/300*stock_concentration)

#get mean concentration per dilution
summarized_intensity_values <- intensity_values %>%
  group_by(Concentration)%>%
  summarize(Intensity=mean(Mean))

#plot trendline
p <- intensity_values %>%
  ggplot(aes(x=Concentration,y=Mean))+geom_point()+  geom_smooth(method=lm) + xlab('Concentration (nM)')
p

#fit data points with linear model extract slope and intercept
m <- lm(Intensity ~Concentration , summarized_intensity_values)

slope <- m$coefficients[2]
intercept <- m$coefficients[1]

#remove too small cells from dataset


#calculate concentration of BRCA2-Halo and number of molecule per nucleus
nucleus_data <- nucleus_data %>%
  mutate(IntegratedIntensity16bit=Intensity_IntegratedIntensity_BRCA2*2^16, Mean_PixInt16bit=IntegratedIntensity16bit/AreaShape_Volume) %>%
  mutate(ConcentrationM=(Mean_PixInt16bit-intercept)/slope*1E-9) %>%
  mutate(Number_of_molecules=ConcentrationM*(AreaShape_Volume*Voxel_volume*1E-15)*Nav)


hist(nucleus_data$ConcentrationM*1E9,xlab="Concentration (nM)")
hist(nucleus_data$Number_of_molecules,main="Number of molecules/nucleus",xlab="Number of molecules/nucleus",breaks=seq(0,30000,2000))

nucleus_data$data <- "WT"
p <- ggplot(nucleus_data,aes(y=ConcentrationM*1E9,x=data))+geom_boxplot()+ 
  geom_quasirandom(width=0.3)+scale_colour_Publication()+scale_fill_Publication()+
  theme_Publication(base_size=12)+ylab("Concentration/nucleus (nM)")
p

p <- ggplot(nucleus_data,aes(y=Number_of_molecules,x=data))+geom_boxplot()+  
  geom_quasirandom(width=0.3)+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=16)+ylab("Number of molecules/nucleus")
p
ggsave(p,filename = file.path(folder,"Nmol_nucleus.png"))

#calculate concentration of BRCA2-Halo and number of molecule per focus
foci_data <- foci_data %>%
  mutate(IntegratedIntensity16bit=Intensity_IntegratedIntensity_BRCA2*2^16, Mean_PixInt16bit=IntegratedIntensity16bit/AreaShape_Volume) %>%
  mutate(ConcentrationM=(Mean_PixInt16bit-intercept)/slope*1E-9) %>%
  mutate(Number_of_molecules=ConcentrationM*(AreaShape_Volume*Voxel_volume*1E-15)*Nav)


hist(foci_data$ConcentrationM*1E9,xlab="Concentration (nM)")
mean(foci_data$AreaShape_Volume)

foci_data$data <- "WT"
ggplot(foci_data,aes(y=Number_of_molecules,x=data))+geom_boxplot(outlier.shape = NA)+ylim(0,75)+
  geom_quasirandom(width=0.3,size=0.2,alpha=0.8)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+
  ylab("Number of molecules/focus")
hist(filter(foci_data,Number_of_molecules<200)$Number_of_molecules,breaks = seq(0,200,10),main="Number of molecules/focus",xlab="Number of molecules/focus")

mean(foci_data$Number_of_molecules)



