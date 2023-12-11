#=====Lecture 8 r code
setwd("C:/Users/timmy/Desktop/FISH_552-3/lecNotes")

#load data
speciesCode <- read.csv("speciesCode.csv")
speciesData <- read.csv("speciesData.csv")
tripData <- read.csv("tripData.csv")
tripData$Date

#tripData examination
head(tripData, n=3)
tripData$Date <- as.Date(tripData$Date)
( min.date <- min(tripData$Date) )      #start date
( max.date <- max(tripData$Date) )      #end date 
( range.date <- range(tripData$Date) )  #min and max dates
difftime(max.date, min.date)
hist(tripData$Date, breaks="months", col="grey") #hist for day data, see ??hist.date

#===longest time difference in data
#order the tripData by the Date column
tripData <- tripData[order(tripData$Date),] #order by date
diff(c(1,2,4,5,6))

#=======In-class exercise 1

max(diff(tripData$Date))
#165 days is the biggest difference

index <- which.max(diff(tripData$Date))
tripData[(index-3):(index+3),] #3 rows above and 3 rows below row 258 with the biggest date difference

#===visualizing trip dates
plot(x=tripData$Date, y=tripData$TotalMinutes/60, 
     type="h", ylim=c(0,5), xaxs="i", yaxs="i",
     xlab="Trip date", ylab="Trip length (hr)")

#===species code data
head(speciesCode, n=4)
head(speciesData, n=4)

#===most frequently caught species
?table
speciesCounts <- table(speciesData$SpeciesCode) #summarizing the number of each species by code
speciesCounts
sort(speciesCounts)
maxCode <- as.numeric(names(speciesCounts[which.max(speciesCounts)]))
maxCode
max.spp <- speciesCode[speciesCode$SpeciesCode==maxCode,]
max.spp


#=====bocaccio dataset
grep("a", c("a","b","a","c","a","d"))
grep("a", c("1a","b","2a","c","3a","d")) 
bocaccioRows <- grep("Bocaccio",speciesCode$Common)
bocaccioRows
speciesCode[bocaccioRows,]

bocaccioCode <- speciesCode[bocaccioRows, "SpeciesCode"]
names(speciesData)
bocaccioData <- subset(speciesData, 
      SpeciesCode==bocaccioCode)
head(bocaccioData)

#merge with the trip data 
bocTrip <- merge(bocaccioData, tripData[,-1],
                 by.x="TripNum", by.y="SimplifiedTripNum")
head(bocTrip, n=4)

#1. Create a data frame by subsetting speciesData to include all rockfish species (scientific name Sebastes) using speciesCode to find the corresponding codes
#2. provide a table() of the fates of rockfish by species code
#3. Calculate the minimum and maximum length recorded for each rockfish species by species code
#. Useful functions include grep(), %in%, table(), names() and tapply(). Also sort(), unique().
#4. Optional (*Advanced): repeat question 2, but obtain a table() of fate vs. common name

#=======In-class exercise 2
#==Q1
rock = speciesCode[275:284,]
rockcode = as.numeric(rock$SpeciesCode)
rockdata = subset(x = speciesData, SpeciesCode %in% rockcode)

#==Q2
table(rockdata$SpeciesCode, rockdata$Fate)

#==Q3
fac = as.factor(unique(rockdata$SpeciesCode))
head(rockdata)
tapply(X = rockdata$Length, as.factor((rockdata$SpeciesCode)), FUN = range, na.rm=TRUE)
tapply(X = rockdata$Length, as.factor((rockdata$SpeciesCode)), FUN = min, na.rm=TRUE)

#==Q4 advanced
r1 = merge(x = rock, y = speciesData[,-1], by.x = "SpeciesCode", by.y = "SpeciesCode")
table(r1$Common, r1$Fate)
