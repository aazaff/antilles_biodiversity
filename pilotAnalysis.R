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
# Steve Holland's optimized SQS function, with error bars added
subsampleEvenness<-function(Abundance,Quota=0.9,Trials=100,IgnoreSingletons=FALSE,ExcludeDominant=FALSE) {
	Abundance<-Abundance[Abundance>0]
	if ((Quota <= 0 || Quota >= 1)) {
		stop('The SQS Quota must be greater than 0.0 and less than 1.0')
		}
 	# compute basic statistics
	Specimens<-sum(Abundance)
	NumTaxa<-length(Abundance)
	Singletons<-sum(Abundance==1)
	Doubletons<-sum(Abundance==2)
	Highest<-max(Abundance)
	MostFrequent<-which(Abundance==Highest)[1]
 	if (ExcludeDominant==FALSE) {
		Highest<-0
		MostFrequent<-0
		}
 	# compute Good's u
	U<-0
	if (ExcludeDominant==TRUE) {
		U<-1-Singletons/(Specimens-Highest)
		} 
	else {
		U<-1-Singletons/Specimens
		}
 	if (U==0) {
		stop('Coverage is zero because all taxa are singletons')
		}
 	# re-compute taxon frequencies for SQS
	FrequencyInitial<-Abundance-(Singletons+Doubletons/2)/NumTaxa
	Frequency<-FrequencyInitial/(Specimens-Highest)
 	# return if the quorum target is higher than estimated coverage
	if ((Quota>sum(Frequency)) || (Quota >= sum(Abundance))) {
		stop('SQS Quota is too large, relative to the estimated coverage')
		}
 	# create a vector, length equal to total number of specimens,
	# each value is the index of that species in the Abundance array
	IDS<-unlist(mapply(rep,1:NumTaxa,Abundance))
 	# subsampling trial loop
	Richness<-rep(0,Trials) # subsampled taxon richness
	for (Trial in 1:Trials) {
		Pool<-IDS # pool from which specimens will be sampled
		SpecimensRemaining<- length(Pool) # number of specimens remaining to be sampled
		Seen<-rep(0,NumTaxa) # keeps track of whether taxa have been sampled
		SubsampledFrequency<-rep(0,NumTaxa) # subsampled frequencies of the taxa
		Coverage<-0
 		while (Coverage<Quota) {
			# draw a specimen
			DrawnSpecimen<-sample(1:SpecimensRemaining,size=1)
			DrawnTaxon<-Pool[DrawnSpecimen]
 			# increment frequency for this taxon
			SubsampledFrequency[DrawnTaxon]<-SubsampledFrequency[DrawnTaxon]+1
 			# if taxon has not yet been found, increment the coverage
			if (Seen[DrawnTaxon]==0) {
				if (DrawnTaxon!=MostFrequent&&(IgnoreSingletons==0||Abundance[DrawnTaxon]>1)) {
					Coverage<-Coverage+Frequency[DrawnTaxon]
					}
				Seen[DrawnTaxon]<-1
 				# increment the richness if the Quota hasn't been exceeded,
				# and randomly throw back some draws that put the coverage over Quota
				if (Coverage<Quota || runif(1)<=Frequency[DrawnTaxon]) {
					Richness[Trial]<-Richness[Trial]+1
					} 
				else {
					SubsampledFrequency[DrawnTaxon]<-SubsampledFrequency[DrawnTaxon]-1
					}
				}
 			# decrease pool of specimens not yet drawn
			Pool[DrawnSpecimen]<-Pool[SpecimensRemaining]
			SpecimensRemaining<-SpecimensRemaining-1
			}
		}
 	# compute subsampled richness
	S2<-Richness[Richness>0]
	SubsampledRichness<-exp(mean(log(S2)))*length(S2)/length(Richness)
	return(c(round(SubsampledRichness,1),quantile(S2,0.975),quantile(S2,0.025)))
	}

############################################# LDG SCRIPT, CHIPS #############################################
# Break data into 1 degree bins (this could be done with just floor(), but round_any is more general)
CanonicalOBIS[,"lat_bin"]<-plyr::round_any(CanonicalOBIS[,"latitude"],1,f=floor)
AbundanceOBIS<-velociraptr::abundanceMatrix(CanonicalOBIS,"genus_name","lat_bin")

# Remove degrees of latitude that are too poorly sampled (Below 1000 individuals)
SizeOBIS<-apply(AbundanceOBIS,2,sum)
AbundanceOBIS<-AbundanceOBIS[,-which(SizeOBIS<1000)]

# Calculate the SQS
BiodiversityOBIS<-pbapply(AbundanceOBIS,2,subsampleEvenness,0.75)

# Make the basic plot
quartz(width=9,height=6)
par(oma=c(1.5,0.5,0.5,0),mar=c(3,3,2,0.5),mgp=c(1.5,0.5,0))
plot(y=BiodiversityOBIS[1,],x=as.numeric(colnames(BiodiversityOBIS)),pch=16,cex=1.5,las=1,xlab="latitude",ylab="standardized richness",xaxs="i",yaxs="i",xlim=c(-90,90),ylim=c(0,525))
arrows(x0=as.numeric(colnames(BiodiversityOBIS)),y0=BiodiversityOBIS[2,],x1=as.numeric(colnames(BiodiversityOBIS)),y1=BiodiversityOBIS[3,],code=0)
abline(v=17); abline(v=19);

# Make a density plot
quartz(width=9,height=6)
par(oma=c(1.5,0.5,0.5,0),mar=c(3,3,2,0.5),mgp=c(1.5,0.5,0))
DistributionOBIS<-rep(as.numeric(colnames(BiodiversityOBIS)),times=ceiling(BiodiversityOBIS[1,]))
plot(density(DistributionOBIS),type="n",xlab="latitude",xlim=c(-90,90),ylab="kernel density of standardized richness",xaxs="i",yaxs="i")
polygon(density(DistributionOBIS),col="cyan",lty=0)
abline(v=17); abline(v=19)

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
                  
# Plot the overall beta diversity consistency
plot(y=TaxonBeta[,"B"],x=TaxonBeta[,"Live"],pch=16,cex=1.5,col="#57afe1",xaxs="i",yaxs="i",xlim=c(0,1),ylim=c(0,1),xlab="life assemblage endemism",ylab="death assemblage endemism")
points(y=TaxonBeta[,"C"],x=TaxonBeta[,"Live"],pch=16,cex=1.5,col="#47a870")
points(y=TaxonBeta[,"D"],x=TaxonBeta[,"Live"],pch=16,cex=1.5,col="#db4544")
		  
# Plot the overall

                  
                  
