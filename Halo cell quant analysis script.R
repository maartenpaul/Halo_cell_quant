#load packages
library(tidyverse)
library(ggbeeswarm)
library(gridExtra)

#modified from https://github.com/koundy/ggplot_theme_Publication
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
            legend.margin = unit(0.1, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(1,1,1,1),"mm"),
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
folder <- "" #input folder
intensity_values <- read_delim(file.path(folder,".txt")) #calibration text file
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
p_trendline <- intensity_values %>%
  ggplot(aes(x=Concentration,y=Mean))+geom_point(size=0.8)+  geom_smooth(method=lm) + xlab('Concentration (nM)')+ylab("Mean intensity")+
  ylim(0,1500)+xlim(0,15)+scale_colour_Publication()+scale_fill_Publication()+
    theme_Publication(base_size=16)
p_trendline
ggsave(p_trendline,filename = file.path(folder,"trendline.png"),width = 10,height=15,units = "cm")
ggsave(p_trendline,filename = file.path(folder,"trendline.pdf"),width = 10,height=15,units = "cm")


#fit data points with linear model extract slope and intercept
m <- lm(Intensity ~Concentration , summarized_intensity_values)

slope <- m$coefficients[2]
intercept <- m$coefficients[1]

#calculate concentration of BRCA2-Halo and number of molecule per nucleus
nucleus_data <- nucleus_data %>%
  mutate(IntegratedIntensity16bit=Intensity_IntegratedIntensity_BRCA2*2^16, Mean_PixInt16bit=IntegratedIntensity16bit/AreaShape_Volume) %>%
  mutate(ConcentrationM=(Mean_PixInt16bit-intercept)/slope*1E-9) %>%
  mutate(Number_of_molecules=ConcentrationM*(AreaShape_Volume*Voxel_volume*1E-15)*Nav)

hist(nucleus_data$ConcentrationM*1E9,xlab="Concentration (nM)")
hist(nucleus_data$Number_of_molecules,main="Number of molecules/nucleus",xlab="Number of molecules/nucleus",breaks=seq(0,30000,2000))

nucleus_data$data <- "WT"
p_nucleus_conc <- ggplot(nucleus_data,aes(y=ConcentrationM*1E9,x=data))+geom_boxplot()+ 
  geom_quasirandom(width=0.3)+scale_colour_Publication()+scale_fill_Publication()+
  theme_Publication(base_size=20)+ylab("Concentration/nucleus (nM)")+ylim(0,25)
p_nucleus_conc
ggsave(p_nucleus_conc,filename = file.path(folder,"Concentration_nucleus.png"),width = 10,height=15,units = "cm")
ggsave(p_nucleus_conc,filename = file.path(folder,"Concentration_nucleus.pdf"),width = 10,height=15,units = "cm")

p_nucleus_n <- ggplot(nucleus_data,aes(y=Number_of_molecules,x=data))+geom_boxplot()+  
  geom_quasirandom(width=0.3)+scale_colour_Publication()+
  scale_fill_Publication()+theme_Publication(base_size=20)+ylab("Number of molecules/nucleus")
p_nucleus_n
ggsave(p_nucleus_n,filename = file.path(folder,"Nmol_nucleus.pdf"),width = 10,height=15,units = "cm")
ggsave(p_nucleus_n,filename = file.path(folder,"Nmol_nucleus.png"),width = 10,height=15,units = "cm")

#calculate concentration of BRCA2-Halo and number of molecule per focus
foci_data$data <- "WT"
foci_data <- foci_data %>%
  mutate(IntegratedIntensity16bit=Intensity_IntegratedIntensity_BRCA2*2^16, Mean_PixInt16bit=IntegratedIntensity16bit/AreaShape_Volume) %>%
  mutate(ConcentrationM=(Mean_PixInt16bit-intercept)/slope*1E-9) %>%
  mutate(Number_of_molecules=ConcentrationM*(AreaShape_Volume*Voxel_volume*1E-15)*Nav)

p_focus_conc <- ggplot(foci_data,aes(y=ConcentrationM*1E9,x=data))+geom_boxplot(outlier.shape = NA)+ylim(0,300)+
  geom_quasirandom(width=0.3,size=0.2,alpha=0.8)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=20)+
  ylab("Concentration/focus (nM)")
p_focus_conc
ggsave(p_focus_conc,filename = file.path(folder,"Concentration_focus.pdf"),width = 10,height=15,units = "cm")
ggsave(p_focus_conc,filename = file.path(folder,"Concentration_focus.png"),width = 10,height=15,units = "cm")
hist(foci_data$ConcentrationM*1E9,xlab="Concentration (nM)")
mean(foci_data$AreaShape_Volume)

p_focus_n <- foci_data %>%
#  dplyr::filter(AreaShape_Volume>10) %>%
    ggplot(aes(y=Number_of_molecules,x=data))+geom_boxplot(outlier.shape = NA)+ylim(0,75)+
  geom_quasirandom(width=0.3,size=0.2,alpha=0.8)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=20)+
  ylab("Number of molecules/focus")
p_focus_n
ggsave(p_focus_n,filename = file.path(folder,"Nmol_focus.pdf"),width = 10,height=15,units = "cm")
ggsave(p_focus_n,filename = file.path(folder,"Nmol_focus.png"),width = 10,height=15,units = "cm")
hist(filter(foci_data,Number_of_molecules<200)$Number_of_molecules,breaks = seq(0,200,10),main="Number of molecules/focus",xlab="Number of molecules/focus")

mean(foci_data$Number_of_molecules)

#combine plots in one overview plot
p_grid <- grid.arrange(p_trendline,p_nucleus_conc,p_nucleus_n,p_focus_conc,p_focus_n,layout_matrix = rbind(c(1, 2, 4),
                                                                                                  c(1, 3, 5)))

