library(ggplot2)
library(functions)

#dir = '../Output_lss'
dir = '../Output_halos'
#dir = '../Output_lss2'


df <- read.csv(paste0(dir,'/kappaAve_z2.297000PS.csv'))

n <- nrow(df)
#df$PS <- df$PS - df$PS[n]
df$llP <- df$l*df$l*df$PS
df <- subset(df,l>90)

df$zs <- '2.297'

#plt <- plt + geom_line(data=df,aes(x=l,y=llP,color='LSS'))

dft <- read.csv(paste0(dir,'/kappaAve_z1.075000PS.csv'))

n <- nrow(dft)
#df$PS <- df$PS - df$PS[n]
dft$llP <- dft$l*dft$l*dft$PS
dft <- subset(dft,l>90)

dft$zs <- '1.075'

df <- rbind(df,dft)
dft <- read.csv(paste0(dir,'/kappaAve_z0.489200PS.csv'))

n <- nrow(dft)
#df$PS <- df$PS - df$PS[n]
dft$llP <- dft$l*dft$l*dft$PS
dft <- subset(dft,l>90)

dft$zs <- '0.489'

df <- rbind(df,dft)

plt <- ggplot(df,aes(x=l,y=llP,colour=zs)) +
  scale_x_log10(limit=c(100,1.0e4)) + scale_y_log10() + 
  ggtitle("Power") +
  xlab("l") + ylab( expression( l^2~P(l) ) ) +
  theme(axis.title.x=element_text(face="italic")) +
  geom_line()

#####################################################
#  BigMultiDark spectra
#####################################################

MDpower <- readVipersBMD('10')

#MDpowert <- readVipersBMD('1')
#MDpower <- rbind(MDpower,MDpowert)

#MDpowert <- readVipersBMD('17')
#MDpower <- rbind(MDpower,MDpowert)

#MDpower$llP <- MDpower$llP/pi/8

plt + geom_line(data=MDpower)

#test <- read.csv("particle_mapPowerSpectrum.dat",sep="")
#test$llP <- test$l*test$l*test$P/pi**2/8
#test$zs <- "1.0"
#test <- subset(test,l>100)
#plt + geom_line(data=test)



