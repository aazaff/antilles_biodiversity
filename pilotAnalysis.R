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
        
# Load the velociraptr package
if (suppressWarnings(require("pbapply"))==FALSE) {
        install.packages("pbapply",repos="http://cran.cnr.berkeley.edu/");
        library("pbapply");
        }

# Establish postgresql connection.
Driver<-dbDriver("PostgreSQL") # Establish database driver
# Your Database
Chips<-dbConnect(Driver, dbname = "chips", host = "localhost", port = 5432, user = "zaffos")
Obis<-dbConnect(Driver, dbname = "obis", host = "localhost", port = 5432, user = "zaffos")

#############################################################################################################
######################################### DATA DOWNLOAD FUNCTIONS, CHIPS ####################################
#############################################################################################################
# No functions at this time
# See the gradientData.R script for how to download the OBIS data

########################################## DOWNLOAD SCRIPTS, CHIPS ##########################################
# Download the data files from postgres
DeathAssemblage2011<-transform(dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.death_assemblage_2011;"),row.names=Sites,Sites=NULL)
DeathAssemblage2012<-transform(dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.death_assemblage_2011;"),row.names=Sites,Sites=NULL)
LiveAssemblage2011<-transform(dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.death_assemblage_2011;"),row.names=Sites,Sites=NULL)
LiveAssemblage2012<-transform(dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.death_assemblage_2011;"),row.names=Sites,Sites=NULL)

# Load the metadata for each locality
Localities<-dbGetQuery(Chips,"SELECT * FROM stcroix_pilot.locality_metadata;")

# Bind the datasets together and round up individuals
DeathAssemblage<-ceiling(rbind(DeathAssemblage2011,DeathAssemblage2012))
LiveAssemblage<-ceiling(rbind(LiveAssemblage2011,LiveAssemblage2012))
# Convert the datasets to presence-absence (This will no longer be necessary starting with velociraptr v1.1)
DeathAssemblage[DeathAssemblage>0]<-1
LiveAssemblage[LiveAssemblage>0]<-1

# Load the OBIS data into R
CanonicalOBIS<-dbGetQuery(Obis,"SELECT * FROM obis_invertebrate_data.obis_data_temp")
# Subset to only the super shellies
# SuperShellies<-c("Bivalvia","Gastropoda","Anthozoa")
# CanonicalOBIS<-subset(CanonicalOBIS,CanonicalOBIS[,"class_name"]%in%SuperShellies==TRUE)
# Extract only rows where taxa are matched to an ecoregion (shallow shelf or island taxa)
# CanonicalOBIS<-subset(CanonicalOBIS,is.na(CanonicalOBIS[,"ecoregion"])!=TRUE)
# Extract only rows where taxa are matched to continental shelf
CanonicalOBIS<-subset(CanonicalOBIS,is.na(CanonicalOBIS[,"plate_id"])!=TRUE)
# Remove taxa without a SST
CanonicalOBIS<-subset(CanonicalOBIS,is.na(CanonicalOBIS[,"temperature"])!=TRUE)
# Remove occurrences without a latitude
CanonicalOBIS<-subset(CanonicalOBIS,is.na(CanonicalOBIS[,"latitude"])!=TRUE)
# Remove bad taxonomic names
CanonicalOBIS<-velociraptr::cleanTaxonomy(CanonicalOBIS,"genus_name")

#############################################################################################################
################################################ LDG PLOT, CHIPS ############################################
#############################################################################################################
# No functions at this time

############################################# LDG SCRIPT, CHIPS #############################################
# Break data into 1 degree bins (this could be done with just floor(), but round_any is more general)
CanonicalOBIS[,"lat_bin"]<-plyr::round_any(CanonicalOBIS[,"latitude"],1,f=floor)
AbundanceOBIS<-velociraptr::abundanceMatrix(CanonicalOBIS,"genus_name","lat_bin")

# Remove degrees of latitude that are too poorly sampled (Below 1000 individuals)
SizeOBIS<-apply(AbundanceOBIS,2,sum)
AbundanceOBIS<-AbundanceOBIS[,-which(SizeOBIS<1000)]

# Calculate the SQS
BiodiversityOBIS<-pbapply(AbundanceOBIS,2,velociraptr::subsampleEvenness,0.75)

# Make the basic plot
plot(y=BiodiversityOBIS,x=as.numeric(names(BiodiversityOBIS)),pch=16,cex=1.5,las=1,xlab="latitude",ylab="standardized richness",xaxs="i",yaxs="i",xlim=c(-90,90),ylim=c(0,525))
abline(v=22)

#############################################################################################################
############################################ ADP FUNCTIONS, CHIPS ###########################################
#############################################################################################################
# No functions at this time

################################################ ADP, CHIPS #################################################
# Separate out the death assemblage samples by stratigraphic code
DeadDivisions<-by(DeathAssemblage,substring(rownames(DeathAssemblage),5,5),function(x) cullMatrix(x,1,1))
# Calculate each taxon's contribution to beta diversity for each stratigraphic layer
DeadBeta<-lapply(DeadDivisions,taxonBeta)

# Merge the data together
TaxonBeta<-transform(merge(DeadBeta[["B"]],DeadBeta[["C"]],by="row.names",all=TRUE),Row.names=NULL,row.names=Row.names)
TaxonBeta<-transform(merge(TaxonBeta,DeadBeta[["D"]],by="row.names",all=TRUE),Row.names=NULL,row.names=Row.names)
                  
# Merge with the live assemblage taxon betas
TaxonBeta<-transform(merge(taxonBeta(cullMatrix(LiveAssemblage,1,1)),TaxonBeta,by="row.names",all=TRUE),Row.names=NULL,row.names=Row.names)
colnames(TaxonBeta)<-c("Live","B","C","D")
                  
# Find the correlation matrix for the beta diversity
cor(TaxonBeta,use="pairwise.complete.obs") 
                  
                  
