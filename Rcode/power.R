df <- read.csv('../Output/kappaAve_z2.297000PS.csv')


df$PS <- df$PS - df$PS[50]
df$llP <- df$l*df$l*df$PS
plot(log10(df$l),log10(df$llP),type='l',col='red')

df <- read.csv('../Output/kappaAve_z0.489200PS.csv')
df$PS <- df$PS - df$PS[50]
df$llP <- df$l*df$l*df$PS

lines(log10(df$l),log10(df$llP),col='blue')


