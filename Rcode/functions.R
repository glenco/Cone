
readVipersBMD <- function(number){
  dir = '/Users/bmetcalf/Projects/'
  
  MDpower <- read.table(paste0(dir,"w4vipers/1/lensing2.0.",number,"/mapPowerSpectrum.dat"),
                      header=FALSE,sep="")
  colnames(MDpower) <- c("l","PS")

  for (i in 2:99) {
    p <- read.table(paste0(dir,"w4vipers/",i,"/lensing2.0.",number,"/mapPowerSpectrum.dat"),
                  header = FALSE,sep="")
  
    MDpower$PS <- MDpower$PS + p$V2
  }
  MDpower$PS <- MDpower$PS/99.
  MDpower$llP <- pi*MDpower$PS*MDpower$l**2
  MDpower <- subset(MDpower,l>100)
  
  return(MDpower)
}