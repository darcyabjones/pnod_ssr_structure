library(ggplot2) # must be version > 2.2.0

# import data
weather <- read_csv("data/weatherdata.csv")
# show first observations
head(weather)
weather$Avearge_Rainfall2 <- weather$Avearge_Rainfall/5



weatherplot <- ggplot(weather,aes(year)) +
      geom_line(aes(y= Min_Temperature, colour = "Min_Temperature"), linetype = "dashed") +
      geom_line(aes(y= Max_Temperature, colour = "Max_Temperature"), linetype = "dashed") + 
      geom_line(aes(y= Avearge_Rainfall/5, colour = "Avearge_Rainfall"),linetype="dashed")+
      geom_smooth(method = "lm", se=FALSE, color="orange", aes(y = Min_Temperature)) +
      geom_smooth(method = "lm", se=FALSE, color="black", aes(y = Max_Temperature)) +
      geom_smooth(method = "lm", se=FALSE, color="red", aes(y = Avearge_Rainfall2))


# now adding the secondary axis, following the example in the help file ?scale_y_continuous
# and, very important, reverting the above transformation
weatherplot <- weatherplot + scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Avearge_Rainfall"))   

# modifying colours and theme options
weatherplot <- weatherplot + scale_colour_manual(values = c("blue", "purple", "dark green"))
weatherplot <- weatherplot + labs(y = "Temperature [°C]", x = "Year", colour = "Parameter")

weatherplot

