library(dataRetrieval) # needed to get peak streamflows to determined active gages
library(kernlab) # for the SVM
library(sp) # Load the sp package for Spatial Object work
library(rgdal) # reading or writing shapefile as needed
library(GISTools) # transparency adder (as written the script does not use add.alpha)
library(feather) # flaky Internet connections can yield problems in repeated pulling
# of USGS peak data, we will write feather files as a type of cache

"insertWaterYear" <- function(x) {
   if(!is.data.frame(x)) {
    stop(paste0("a data.frame as in 'dataRetrieval::readNWISpeak(station, ",
                "convertType=FALSE)' required"))
  }
  x$peak_va  <- as.numeric(x$peak_va)
  x$year_va  <- as.numeric(sapply(x$peak_dt,
                                  function(t) strsplit(t, split = "-")[[1]][1]))
  x$month_va <- as.numeric(sapply(x$peak_dt,
                                  function(t) strsplit(t, split = "-")[[1]][2]))
  x$day_va   <- as.numeric(sapply(x$peak_dt,
                                  function(t) strsplit(t, split = "-")[[1]][3]))
  x$water_yr <- x$year_va; getem <- ! is.na(x$month_va) & x$month_va >= 10
  x$water_yr[getem] <- x$water_yr[getem] + 1
  return(x)
}


"myreadNWISpeak" <- function(sites, feathercache=TRUE) {
   if(length(sites) == 0) return(NULL)
   first <- TRUE
   zz <- NULL
   for(i in 1:length(sites)) {
      message(i,"-", appendLF=FALSE)
      pkfile <- paste0("pkrda/",sites[i],".feather")
      pk <- NULL
      if(feathercache & file.exists(pkfile)) {
        pk <- read_feather(pkfile)
      } else {
        try(pk <- dataRetrieval::readNWISpeak(sites[i], convertType=FALSE), silent=TRUE)
      }
      if(is.null(pk)) {
         message("failed ",sites[i])
         next
      }
      if(length(pk) == 1) {
         message("failed ",sites[i])
         next
      }
      pk <- insertWaterYear(pk)
      if(is.null(pk)) next # could be "peaks" in the database are all NA
      if(feathercache) write_feather(pk, pkfile)
      if(first) { zz <- pk; first <- FALSE; next }
      zz <- rbind(zz, pk)
      #print(head(pk, n=2))
   }
   message("done")
   if(is.null(zz)) return(NULL)
   return(zz)
}


LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG) # lat and longitude (geographic coordinate definition)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA) # Albers Equal Area conic projection definition

ST <- rgdal::readOGR("data/cb_2017_us_state_500k", "cb_2017_us_state_500k")
st <- slot(ST, "polygons")
#for(i in 1:50) {
TX <- sp::SpatialPolygons(list(st[[37]]), proj4string=CRS(proj4string(ST)))
#plot(TX); mtext(i)
#}
TX <- sp::spTransform(TX, ALBEA); rm(ST); rm(st)

WC <- read.table("data/Appendix1_638wtrshdchr.txt",
                 sep=",", header=TRUE, colClasses="character")
QT <- read.table("data/Appendix1_677trimmedQTs.txt",
                 sep=",", header=TRUE, colClasses="character")
WC$CDA <- as.numeric(WC$CDA)
WC$MCS <- as.numeric(WC$MCS)
WC$MAP <- as.numeric(WC$MAP)
QT$trimQ10 <- as.numeric(QT$trimQ10)

DB <- merge(WC, QT, by.x="STATION", all=TRUE)
DB <- DB[complete.cases(DB),]
DB <- DB[DB$STATION != "08098203",] # NWIS site file is lacking this gage
DB <- DB[DB$STATION != "08086300",] # NWIS does not seem to be now (04/19/2019) reporting data
print(length(unique(DB$STATION)))
print(length(DB$STATION))

Q10 <- log10(DB$trimQ10); CDA <- log10(DB$CDA)
MCS <- log10(DB$MCS);     MAP <- log10(DB$MAP)


SF <- dataRetrieval::readNWISsite(DB$STATION)
SF <- SF[SF$agency_cd != "USIBW",] # remove International Boundary Water Commission sites
SF <- sp::SpatialPointsDataFrame(cbind(SF$dec_long_va,SF$dec_lat_va),
                             data=data.frame(STATION=SF$site_no,
                                             NAME=SF$station_nm,
                                             STATE=SF$state_cd,
                                             LAT=SF$dec_lat_va, LON=SF$dec_long_va,
                                             CDA=CDA, MCS=MCS, MAP=MAP, Q10=Q10),
                             proj4string=LATLONG)
SF <- sp::spTransform(SF, ALBEA); XY <- coordinates(SF)
SF$X <- XY[,1]/1000; SF$Y <- XY[,2]/1000; rm(XY)
SF <- SF[SF$STATE == 48,] # isolate only Texas

print(length(SF$STATION)) # [1] 536 streamgages


message("Summary of Overall Network Contributing Drainage Areas in km^2")
print(10^summary(SF$CDA)*2.589988)
message("Summary of Overall Network Main-Channel Slopes, dimensionless")
print(10^summary(SF$MCS))
message("Summary of Overall Network Mean Annual Precipitation in mm")
print(10^summary(SF$MAP)*25.4)

# Begin the simulations
nsim <- 20
first <- TRUE; ix <- 1:length(SF$Q10)
message("Starting simulation chucks totaling ",nsim," and this takes time!")
for(j in 1:nsim) {
  set.seed(j)
  message("simulation chuck: ",j)
  for(i in 1:length(SF$Q10)) {
    if(as.logical(length(grep("00$", i)))) message(i, "-", appendLF=FALSE)
    SVM <- kernlab::ksvm(SF$Q10[-i]~SF$CDA[-i]+SF$MCS[-i]+SF$MAP[-i]+SF$X[-i]+SF$Y[-i])
    tix <- ix[-i]; six <- tix[SVindex(SVM)]
    if(first) { first <- FALSE
      svm <- data.frame(index=six)
    } else {
      svm <- rbind(svm, data.frame(index=six))
    }
  }
  message("done")
}
svm <- aggregate(svm, by=list(svm$index), length)
svm$svm_count <- svm$index; svm$index <- svm$Group.1; svm$Group.1 <- NULL
svm$STATION <- SF$STATION[svm$index]
svm$svm_ratio <- svm$svm_count/(nsim*(length(SF$Q10)-1))

# Now join the ratios into the sitefile (SF) data structure
SF$svm_ratio <- SF$svm_count <- 0
for(site in SF$STATION) {
  n <- length(svm$svm_count[svm$STATION == site])
  message(n,"-", appendLF=FALSE)
  if(n == 0) { message(""); next }
  SF$svm_count[SF$STATION == site] <- svm$svm_count[svm$STATION == site]
  SF$svm_ratio[SF$STATION == site] <- svm$svm_ratio[svm$STATION == site]
}



# In LaTeX source \numgageslessfivepercentSVs is the variable less_than_5
less_than_five <- length(SF$STATION[SF$svm_ratio <= 0.05])
# In LaTeX source \percentgageslessfivepercentSVs is the variable pct_less_than_five
pct_less_than_five <- round(100*less_than_five/length(SF$STATION), digits=0)
message("Percent of streamgages <= 5 percent of time support vectors = ", pct_less_than_five," percent")
message("Number of streamgages <=5 percent of time support vectors = ", less_than_five)


needed_greatly <- as.character(SF$STATION[SF$svm_ratio == 1])
always_hundred_percent <- length(needed_greatly)
# In LaTeX source \numgagesalwaysSVs is the variable always_hundred_percent
message("Number of streamgages 100 percent of time support vectors = ", always_hundred_percent)

txt <- "Color hue is prorated from red to blue\nbased on nonexceedance probability."
pdf("../draftfigures/fig09_svmtexaspp.pdf", useDingbats=TRUE)
  opts <- par(no.readonly = TRUE); par(las=1)
  tmp <- SF[order(SF$svm_ratio),]
  plot(lmomco::pp(tmp$svm_ratio, sort=FALSE),
                  tmp$svm_ratio*100, xlim=c(0,1),
       xlab="NONEXCEEDANCE PROBABILITY", lwd=0.8, type="p",
       ylab="SUPPORT VECTOR PERCENTAGE",
       col=rgb(1-tmp$svm_ratio,0,tmp$svm_ratio))
  text(0.7, 50, txt, cex=0.9)
  par(opts)
dev.off()


message("Now going to the Internet to pull peaks to determine 'active'",
        " or 'discontinued' according to settings of this algorithm")
if(! dir.exists("./pkrda/")) dir.create("./pkrda/")
if(file.exists("./pkrda.zip")) unzip("./pkrda.zip", overwrite=TRUE)
PK <- myreadNWISpeak(SF$STATION)
if(dir.exists("./pkrda/"))   zip("./pkrda.zip", "./pkrda")
if(dir.exists("./pkrda/"))   unlink("./pkrda/",   recursive=TRUE)
if(dir.exists("./__MACOSX")) unlink("./__MACOSX", recursive=TRUE)


last_year <- rep(NA, length(needed_greatly))
i <- 0
for(site in needed_greatly) {
  i <- i + 1
  #pk <- get(site, envir=PKenv)
  pk <- PK[PK$site_no == site,]
  pk <- pk[pk$year_va <= 2017,] # On 10/03/2019 last year was still 2017
  # USGS is behind in updating the peak values database during database
  # transitions. By setting 2017 year, we get reproducibility of results
  # of sources in this code base for purposes of a research paper.
  last_year[i] <- max(pk$year_va)
}

needed_badly <- needed_greatly[last_year < 2000]
discontinued_hundred_percent <- length(needed_badly)
# In LaTeX source \numgagesalwaysSVsdiscontinued is the variable discontinued_hundred_percent
message("Number of streamgages 100 percent of time support vectors that are ",
        "discontinued \\numgagesalwaysSVsdiscontinued= ",discontinued_hundred_percent)

txt <- "Color hue is prorated from\nred to blue based on\nnonexceedance probability."

# Create the map
pdf("../draftfigures/fig10_svmtexasmap.pdf", useDingbats=TRUE)
  opts <- par(no.readonly = TRUE); par(las=1)
  plot(TX, lwd=1.1)
  #plot(SF, lwd=0.5, pch=4, cex=0.8, add=TRUE)
  tmp <- SF[SF$svm_ratio > 0.05 & SF$svm_ratio < 1,]
  plot(tmp, pch=21, col=1, cex=0.6, bg=rgb(1-tmp$svm_ratio,0,tmp$svm_ratio), lwd=0.6, add=TRUE)
  plot(SF[SF$svm_ratio == 1,], pch=24, col=8, bg=4, lwd=0.6, cex=0.8, add=TRUE)
  plot(SF[SF$svm_ratio <= 0.05,], pch=25, col=8, bg=2, lwd=0.6, cex=0.8, add=TRUE)
  for(site in needed_badly) {
    plot(SF[SF$STATION == site,], pch=1, col=4, lwd=0.8, cex=1.2, add=TRUE)
  }
  legend(-340000, 1500000,
         c("ASV sites (SV 100 percent of the time)",
           "NSV sites (SV <=5 percent of the time)",
           "Neither ASV or NSV sites with color blue to red color ramp",
           "Indicator of ASV not operated since at least 2000"),
          bty="n", cex=0.6, pt.bg=c(4,2,rgb(0.5,0,0.5),4),
          pt.lwd=c(0.6,0.6,0.6,0.8), pch=c(24,25,23,1), col=c(8,8,8,4))
  text(-280000, 1330000, "SV, support vector", cex=0.6, pos=4)
  text(-280000, 1300000, "ASV, always a support vector", cex=0.6, pos=4)
  text(-280000, 1270000, "NSV, nonsupport vector", cex=0.6, pos=4)
  north.arrow(-810000,1170000, 15000, col="black", cex=0.6)
  par(cex=0.7)
  map.scale(-715000, 570000, 1000*400, "kilometers", 4, 100)
  text(-715000, 450000, paste0("Base map from U.S. Census Bureau digital sources, 1:500k\n",
                               "Albers equal area projection"),
       cex=0.7)
  text(70000, 450000, txt, cex=0.85)
  par(cex=1)
  par(opts)
dev.off()



message("All \\numgages=",length(SF$svm_ratio),"\n",
        "Num of needed greatly \\numgagesalwaysSVs= ",length(needed_greatly),"\n",
        "Num of needed badly (discontinued sites) \\numgagesalwaysSVsdiscontinued = ",length(needed_badly))

JK <- SF[SF$STATION == needed_badly[1],]
for(i in 2:length(needed_badly)) {
  JK <- rbind(JK,SF[SF$STATION == needed_badly[i],] )
}
message("Summary of ASV Contributing Drainage Areas in km^2")
print(10^summary(JK$CDA)*2.589988)
message("Summary of ASV Main-Channel Slopes, dimensionless")
print(10^summary(JK$MCS))
message("Summary of ASV Mean Annual Precipitation in mm")
print(10^summary(JK$MAP)*25.4)

message("Summary of Overall Network Contributing Drainage Areas in km^2")
print(10^summary(SF$CDA)*2.589988)
message("Summary of ASV Contributing Drainage Areas in km^2")
print(10^summary(JK$CDA)*2.589988)

