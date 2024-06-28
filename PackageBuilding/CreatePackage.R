usethis::create_package('C:/Users/r816g589/Documents/GitHub/micRoclean')

# add .Rd documentation files to man folder
devtools::document('C:/Users/r816g589/Documents/GitHub/micRoclean')

devtools::build_manual('C:/Users/r816g589/Documents/GitHub/micRoclean')

devtools::build('C:/Users/r816g589/Documents/GitHub/micRoclean')

devtools::install_local('C:/Users/r816g589/Documents/GitHub/micRoclean')

# check local installation
devtools::install_local('micRoclean')

# check github installation
devtools::install_github('rachelgriffard/micRoclean')

# remove copy
remove.packages('micRoclean')