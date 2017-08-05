library(ggplot2)
library(functions)

#dirs <- c('../Output_lss','../Output_halos')
#dir = '../Output_lss2'

dir <- "../Output/"
#zs <- "1.075000"
zs <- "0.489200"
labels = c('PowTOTAL','PowHALO','PowLSS')

df <- data.frame()
#norms <- c(1.0,2.5e4)
norms <- c(1.0,1.0,1.0e-0)
i=1
print(dir)
for(label in labels){
  print(label)
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

  file = paste0(dir,'kappaAve_z',zs,label,'.csv')
  print(file)
  dft <- read.csv(file)

  dft$llP <- dft$l*dft$l*(dft$PS - dft$PS[nrow(dft)])/norms[i]
  i = i + 1
  dft <- subset(dft,l>90)

  dft$zs <- paste(zs, label)

  df <- rbind(df,dft)
}

plt <- ggplot(df,aes(x=l,y=llP,colour=zs)) +
  scale_x_log10(limit=c(100,1.0e5)) + scale_y_log10() + 
  ggtitle("Power") +
  xlab("l") + ylab( expression( l^2~P(l) ) ) +
  theme(axis.title.x=element_text(face="italic")) +
  geom_line()


#####################################################
#  BigMultiDark spectra
#####################################################

plt

MDpower <- NULL

MDpowert <- readVipersBMD('18')
MDpower <- rbind(MDpower,MDpowert)

MDpower$llP <- MDpower$llP/pi/8

MDpower$zs <- paste(MDpower$zs,'BigMD')

plt + geom_line(data=MDpower)

#labels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24')
#MDpower <- NULL
#for(i in labels ){
#  print(i)
#  MDpowert <- readVipersBMD(i)
#  MDpower <- rbind(MDpower,MDpowert)
#}

#write.csv(MDpower,file='BigMD.csv')
