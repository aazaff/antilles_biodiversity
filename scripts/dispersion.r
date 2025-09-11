# This is a script of pilot analyses of the St Croix data meant to create a figure that will
# accomplish the following: 1) illustrate the concept of ordination, 2) illustrate some of our community taxa
# and 3) illustrate the concept of "dispersion"

# Some notes on the code, I have adopted a number of more pythonic conventions
# 1. snake_case for arguments and variables instead fo the PascalCase generally preferred by the R Community
# 2. Using the = operator instead of the <- assigment operator unless if directionality is required
# 3. Using []-notation instead of the $ operator, as [] will work with atomic vectors
# 4. Calling variables by name instead of by position in most cases
# 5. I retain the R convention of using camelCase for custom functions
# 6. This was coded on a mac, which retains some differences in system interaction compared to windows, most
# famously is the use of quartz() vs windows() for opening new windows, and how file paths are named.

######################################### Load Required Libraries ###########################################
# Load or install the velociraptr package; needed to shape the data
if (suppressWarnings(require("velociraptr"))==FALSE) {
    install.packages("velociraptr",repos="http://cran.cnr.berkeley.edu/");
    library("velociraptr");
    }

# Load or install the vegan package; needed for ordination
if (suppressWarnings(require("vegan"))==FALSE) {
    install.packages("vegan",repos="http://cran.cnr.berkeley.edu/");
    library("vegan");
    }

#############################################################################################################
###################################### FOSSIL DATA FUNCTIONS, ST CROIX ######################################
#############################################################################################################
# No functions at this time.

############################################## Load St Croix Datasets  ######################################
# Loaded from a shared dropbox folder with Kelsey Arkle and Andrew Zaffos
raw_data = read.csv("~/Dropbox/Augustana RUI Coaching Project/Data (for Andrew)/St_Croix_Dead_Raw.csv",row.names=1)

# Quickly view the data for consistency and get its dimensions
dim(raw_data)
head(raw_data)

# Remove depauperate samples (<2 species) and singleton taxa (<2 samples)
culled_data = velociraptr::cullMatrix(raw_data,2,2)

#############################################################################################################
####################################### ORDINATION FUNCTIONS, ST CROIX ######################################
#############################################################################################################
# Perform an affine rotation so that the slope of the inferred gradient is zero
# This is useful for maximizing correspondence of the ecological signal with gradient axis 1
# You can determined the slope of the inferred gradient as slope = sd(axis2_centorid)/sd(axis1_centroid)
# From Holland and Zaffos (2011)
rotateMatrix <- function(community, theta) {
	# rotates a 2-dimensional matrix clockwise by theta radians
	trigonometry = c(cos(theta), sin(theta), -sin(theta), cos(theta))
	dimensions = c(2,2) # this is hard-coded, and I don't even remember what it does
	output <-array(trigonometry, dimensions)
	output <-community %*% output
	return(output)
	}

######################################## ORDINATION SCRIPT, ST CROIX ########################################
# Standardize the abundance of taxa within samples as proporotional abundances
stand_data = vegan::decostand(culled_data,"total")

# Generate the initial DCA plot
# Note that I do not down-weight rare taxa - i.e., iweigh=0. I don't see a theoretical reason to do so in this case.
dca_results = vegan::decorana(stand_data)

# plot the initial results
plot(dca_results) # based on qualitative analysis the "true" axis 1 gradient seems to be depicted at a 45° angle. (i.e., slope = 1)

# extract the relevant scores for axes 1 and 2
site_scores = scores(dca_results,display="sites",choices=c(1,2))
species_scores = scores(dca_results,display="species",choices=c(1,2))

# Rotate the species scores accordingly
rotated_species = rotateMatrix(species_scores,theta=atan(1))

# Plot the rotated species scores
plot(rotated_species,col="darkgrey",las=1,xlab="Gradient 1",ylab="Gradient 2",pch=16)

#############################################################################################################
####################################### DISPERSION FUNCTIONS, ST CROIX ######################################
#############################################################################################################
# This function iteratively draws a light dashed line using arrows() from centroid to each community member
# and also outputs the average distance. This required plot.new() to have been called already
centroidArrows <- function(community,color="lightgrey") {
    centroid_x = mean(community[,1])
    centroid_y = mean(community[,2])
    for (i in 1:nrow(community)) {
        arrows(centroid_x,centroid_y,community[i,1],community[i,2],lty=3,col=color,code=0,lwd=1.25)
        }
    }

######################################## DISPERSION SCRIPT, ST CROIX ########################################
# Select a "community" of interest. I arbitraritly chose "dump site" for this example, sample code DS
# You can substitute a different locality by just changing the subset string - e.g., "DS" to "SC" for smuggler's cove
dump_site = subset(raw_data,substring(rownames(raw_data),1,2)=="DS")
# I included an optional, commented out approach if you just want to do the "B"  samples of a community
# dump_site = subset(raw_data,substring(rownames(raw_data),1,2)=="DS" & substring(rownames(raw_data),5,5)=="B")

# Again, cull the community of irrelevant samples and species
culled_dump = velociraptr::cullMatrix(dump_site,2,2)

# Extract the species scores for dump_site taxa
dump_community = subset(rotated_species,rownames(rotated_species)%in%colnames(culled_dump)==TRUE)

# (re)add the dump site taxa to the plot
points(dump_community,pch=16,col="black")
# You can substitute the names for the points using text()
# text(dump_community,label=rownames(dump_community))

# calculate the centroid of the community and its dispersion (arithmetic mean distance of species to centroid)
# I wonder if a case could be made that this should be calcualted using the geometric mean?
# Also, as a note, I'm using the most basic method of multivariate dispersion as presented by Anderson et al. 2006,
# NOT the weighted-by-abundance variants described by Laliberte and Legendre 2010. I see no reason for this since the
# abundances have ALREADY been used as part of our ordination process.
centroid = c(mean(dump_community[,1]),mean(dump_community[,2]))
dispersion = mean(apply(dump_community,1,function(x) dist(rbind(x,centroid)))) # this could be done more efficiently if I understood slicing the output of dist better, but I don't so whatever.

# draw lines from the centroid to the community members. I chose lightgrey dashed lines. The dashes are hard-coded
# but the color can be changed using the color argument
centroidArrows(dump_community,color="lightgrey")

# Add the centroid of the dumpsite community
points(x=mean(dump_community[,1]),y=mean(dump_community[,2]),col="red",pch=17,cex=1.25)