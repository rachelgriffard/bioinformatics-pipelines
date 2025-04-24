usethis::create_package('') # path to package

# add .Rd documentation files to man folder
devtools::document('C:/Users/path/to/package/pkgname') # path to package

devtools::build_manual('C:/Users/path/to/package/pkgname') # path to package

devtools::build('C:/Users/path/to/package/pkgname') # path to package

devtools::install_local('C:/Users/path/to/package/pkgname') # path to package

# check local installation
devtools::install_local('pkgname') # package name and tar file

# check github installation
devtools::install_github('user/pkgname') # user/repository

# remove copy
remove.packages('pkgname') # delete