setwd("C:/Users/timmy/Desktop/FISH_552-3/lecNotes")
bats <- read.csv("bats.csv")
head(bats)

bat_analysis <- function(common_species_name, ...) {
  x <- subset(bats, common_name == common_species_name)
  bat_lm <- lm(x$wingspan ~ x$weight)
  plot(x = x$wingspan, y = x$weight, xlab = "wingspan", ylab = "weight", 
       main = common_species_name, ...)
  abline(lm(x$weight ~ x$wingspan))
  p_val <- summary(bat_lm)$coefficients[2,4] #p-val
  pp_val <- paste("p value =", signif(p_val, digits = 3))
  bat_list <- list(bat_lm , x, p_val)
  return(bat_list)
}

fit_bats <- function(common_name) {
PB <- bat_analysis(common_name)
females <- subset(PB[[2]], sex == "F")
fit_females <- subset(females, weight > mean(females$weight, na.rm = TRUE) & wingspan > mean(females$wingspan, na.rm = TRUE))
print(fit_females)
}

for (name in unique(bats$common_name)) {
  if (bat_analysis(name)[[3]] < 0.05) {
    print(paste("the weight and wingspan of", name, "has a significant linear relationship"))
    print("lsited below are females above average in weight and wingspan")
    fit_bats(common_name = name)
  } else {
    print(paste(name, "does not fit the requirement"))
  }
}

#Create a new data frame that contains only bats of the genus Myotis
#Count the number of individuals that we have of each species
#Create a new field that contains the weight to wingspan ratio
#Calculate the mean weight to wingspan ratio for each species. Which species is the chunkiest?
#Order your Myotis data so that the chunkiest bats are at the top. Which bat is the chunkiest bat?

substring(bats$scientific_name, 1, 6)

new <- subset(bats, scientific_name == "Myotis californicus" | "Myotis yumanensis" | "Myotis lucifugus")
