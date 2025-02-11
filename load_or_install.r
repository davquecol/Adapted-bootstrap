#load_or_install function
#Version: 2.0
#Last_update:2024/12/24
#Programmer: David_P_Quevedo

load_or_install <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, repos = "https://cloud.r-project.org")  # Explicitly set CRAN mirror
    library(package, character.only = TRUE)
  }
}

