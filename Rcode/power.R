library(ggplot2)
library(functions)

dirs <- c('../Output_lss','../Output_halos')
#dir = '../Output_lss2'


df <- data.frame()

for(dir in dirs){
  print(dir)
  
#  dft <- read.csv(paste0(dir,'/kappaAve_z2.297000PS.csv'))

#  dft$llP <- dft$l*dft$l*dft$PS
#  dft <- subset(dft,l>90)

#  dft$zs <- paste('2.297',dir)

#  df <- rbind(df,dft)
#  dft <- read.csv(paste0(dir,'/kappaAve_z1.075000PS.csv'))

#  dft$llP <- dft$l*dft$l*dft$PS
#  dft <- subset(dft,l>90)

#  dft$zs <- paste('1.075',dir)

#  df <- rbind(df,dft)

  dft <- read.csv(paste0(dir,'/kappaAve_z0.489200PS.csv'))

  dft$llP <- dft$l*dft$l*dft$PS
  dft <- subset(dft,l>90)

  dft$zs <- paste('0.489', dir)

  df <- rbind(df,dft)
}

plt <- ggplot(df,aes(x=l,y=llP,colour=zs)) +
  scale_x_log10(limit=c(100,1.0e4)) + scale_y_log10() + 
  ggtitle("Power") +
  xlab("l") + ylab( expression( l^2~P(l) ) ) +
  theme(axis.title.x=element_text(face="italic")) +
  geom_line()


#####################################################
#  BigMultiDark spectra
#####################################################

MDpower <- NULL

#MDpower <- readVipersBMD('10')

#MDpowert <- readVipersBMD('1')
#MDpower <- rbind(MDpower,MDpowert)

MDpowert <- readVipersBMD('17')
MDpower <- rbind(MDpower,MDpowert)

MDpower$llP <- MDpower$llP/pi/8

MDpower$zs <- paste(MDpower$zs,'BigMD')

plt + geom_line(data=MDpower)




