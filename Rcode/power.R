library(ggplot2)
#library(functions)

plt <- ggplot() +
  scale_x_log10(limit=c(100,1.0e5)) + scale_y_log10() + ggtitle("Power") +
  xlab("l") + ylab( expression( l^2~P(l) ) ) +
  theme(axis.title.x=element_text(face="italic"))


df <- read.csv('../Output/kappaAve_z2.297000PS.csv')

n <- nrow(df)
df$PS <- df$PS - df$PS[n]
df$llP <- df$l*df$l*df$PS
plt <- plt + geom_line(data=df,aes(x=l,y=llP))

df <- read.csv('../Output/kappaAve_z0.489200PS.csv')
n <- nrow(df)
df$PS <- df$PS - df$PS[n]
df$llP <- df$l*df$l*df$PS

plt <- plt + geom_line(data=df,aes(x=l,y=llP,color='red'))

MDpower <- readVipersBMD('10')
plt <- plt + geom_line(data=MDpower
                       ,aes(x=l,y=llPS,color='particles'))

#MDpower <- readVipersBMD('1')
#plt <- plt + geom_line(data=MDpower
#                       ,aes(x=l,y=llPS,color='particles'))

#MDpower <- readVipersBMD('17')
#plt <- plt + geom_line(data=MDpower
#                       ,aes(x=l,y=llPS,color='particles'))

plt
