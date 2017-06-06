
readVipersBMD <- function(number){
  dir = '/Users/bmetcalf/Projects/MultiDPowerSpectrumData/'
  
  z  <- c('2.297','2.119','1.955','1.802','1.66','1.527','1.403','1.287','1.178','1.075','0.9774','0.8854','0.7982','0.7154','0.6365','0.5612','0.4892','0.4201','0.3538','0.2899','0.2282','0.1686','0.1108','0.05465')
  
  MDpower <- read.table(paste0(dir,"1/lensing2.0.",number,"/mapPowerSpectrum.dat"),
                      header=FALSE,sep="")
  colnames(MDpower) <- c("l","PS")

  for (i in 2:99) {
    p <- read.table(paste0(dir,i,"/lensing2.0.",number,"/mapPowerSpectrum.dat"),
                  header = FALSE,sep="")
  
    MDpower$PS <- MDpower$PS + p$V2
  }
  MDpower$PS <- MDpower$PS/99.
  MDpower$llP <- MDpower$PS*MDpower$l**2
  MDpower$zs <- z[as.numeric(number)]
  MDpower <- subset(MDpower,l>100)
  
  return(MDpower)
}