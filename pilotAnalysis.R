# Custom functions are camelCase. DataObjects and Arguments are PascalCase
# Dependency functions are not embedded in master functions
# []-notation is used wherever possible, and $-notation is avoided.

######################################### Load Required Libraries ###########################################
# Allows connection to postgres
if (suppressWarnings(require("RPostgreSQL"))==FALSE) {
        install.packages("RPostgreSQL",repos="http://cran.cnr.berkeley.edu/");
        library("RPostgreSQL");
        }

# Load the velociraptr package
if (suppressWarnings(require("velociraptr"))==FALSE) {
        install.packages("velociraptr",repos="http://cran.cnr.berkeley.edu/");
        library("velociraptr");
        }
        
# Establish postgresql connection.
Driver<-dbDriver("PostgreSQL") # Establish database driver
# Your Database
Chips<-dbConnect(Driver, dbname = "chips", host = "localhost", port = 5432, user = "zaffos")

#############################################################################################################
######################################### DATA DOWNLOAD FUNCTIONS, CHIPS ####################################
#############################################################################################################
# No functions at this time

########################################## DOWNLOAD SCRIPTS, CHIPS ##########################################
# Download the data files from postgres
DeathAssemblage2011<-transform(dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.death_assemblage_2011;"),row.names=Sites,Sites=NULL)
DeathAssemblage2012<-transform(dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.death_assemblage_2011;"),row.names=Sites,Sites=NULL)
LiveAssemblage2011<-transform(dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.death_assemblage_2011;"),row.names=Sites,Sites=NULL)
LiveAssemblage2012<-transform(dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.death_assemblage_2011;"),row.names=Sites,Sites=NULL)

# Bind the datasets together and round up individuals
DeathAssemblage<-ceiling(rbind(DeathAssemblage2011,DeathAssemblage2012))
LiveAssemblage<-ceiling(rbind(LiveAssemblage2011,LiveAssemblage2012))

