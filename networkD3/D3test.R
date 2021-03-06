BiocManager::install("networkD3")
library("networkD3")

# Create fake data
src <- c("A", "A", "A", "A",
         "B", "B", "C", "C", "D")
target <- c("B", "C", "D", "J",
            "E", "F", "G", "H", "I")
networkData <- data.frame(src, target)

# Plot
library(magrittr)

simpleNetwork(networkData,fontSize = 15,nodeColour = "#EA6A47")