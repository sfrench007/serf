# Where the SDF files are in relation to the working directory.  By default here, there is a directory just called 'SDF'
sdfDirectory <- "SDF/"  

# Create a vector of the SDF files to loop through.  If they have capitaliation in the filenames, that needs to be specified below, or the 'pattern' option removed
sdfList <- list.files(sdfDirectory,pattern=".sdf",full.names=TRUE,recursive=FALSE)
# Get a list of the non-full-filename files
sdfNames <- list.files(sdfDirectory,pattern=".sdf",full.names=FALSE,recursive=FALSE)
# Remove the filename to create the .txt file outputs
outputList <- paste(substr(sdfNames,1,nchar(sdfNames)-3),"txt",sep="")

# The chemical properties from cxcalc
props <- c("atomcount",
"exactmass",
"massspectrum",
"atomicpolarizability --pH 7",
"atompol --pH 7",
"axxpol --pH 7",
"ayypol --pH 7",
"azzpol --pH 7",
"charge --pH 7",
"formalcharge --pH 7",
"ioncharge --pH 7",
"oen --pH 7",
"polarizability --pH 7",
"tholepolarizability --pH 7",
"aliphaticatom",
"aliphaticatomcount",
"aliphaticbondcount",
"aliphaticringcount",
"aliphaticringcountofsize",
"asa",
"asymmetricatom",
"asymmetricatomcount",
"balabanindex",
"carboaliphaticringcount",
"carboringcount",
"chainatom",
"chiralcenter",
"chiralcentercount",
"cyclomaticnumber",
"distancedegree",
"eccentricity",
"fsp3",
"fusedaliphaticringcount",
"fusedaromaticringcount",
"hararyindex",
"heteroaliphaticringcount",
"heteroaromaticringcount",
"hyperwienerindex",
"largestatomringsize",
"molecularsurfacearea --pH 7",
"msa --pH 7",
"plattindex",
"polarsurfacearea --pH 7",
"psa --pH 7",
"randicindex",
"rotatablebondcount",
"stereodoublebondcount",
"stericeffectindex",
"szegedindex",
"vdwsa --pH 7",
"wateraccessiblesurfacearea --pH 7",
"wienerindex",
"wienerpolarity",
"doublebondstereoisomercount",
"tautomercount --pH 7",
"tetrahedralstereoisomercount",
"logd --pH 7",
"logp --pH 7",
"averagemicrospeciescharge --pH 7",
"chargedistribution --pH 7",
"isoelectricpoint",
"pka",
"acc --pH 7",
"acceptor --pH 7",
"acceptormultiplicity --pH 7",
"chargedensity --pH 7",
"don --pH 7",
"donor --pH 7",
"donorcount --pH 7",
"donormultiplicity --pH 7",
"electrondensity --pH 7",
"electrophiliclocalizationenergy --pH 7",
"hbda --pH 7",
"hmochargedensity --pH 7",
"hmoelectrondensity --pH 7",
"hmoelectrophilicityorder --pH 7",
"hmoelectrophiliclocalizationenergy --pH 7",
"hmohuckeleigenvalue --pH 7",
"hmohuckelorbitals --pH 7",
"hmolocalizationenergy --pH 7",
"hmonucleophilicityorder --pH 7",
"hmopienergy --pH 7",
"nucleophiliclocalizationenergy --pH 7",
"refractivity",
"resonantcount"
)

# Bind the property arguments together into a single line
props2 <- ""; for(i in 1:length(props)) props2 <- paste(props2,props[i],sep=" ")

# Create a progress bar for the # of SDF files, and call cxcalc to calculate the properties.  Automatically saved as a .txt file with the same filename
pb <- txtProgressBar(min=1,max=length(sdfList),style=3)
for(i in 1:length(sdfList)) {
	print(paste("Analyzing ",sdfNames[i],sep=""))
	sdfFile <- sdfList[i]
	outputFile <- outputList[i]
	setTxtProgressBar(pb,i)
	
	# The system command to call cxcalc.  Chemaxon cxcalc must be installed and *registered* in order for the command to work.  Users will need to specify
	# the full path of cxcalc here if it is not installed globally
	system(paste("cxcalc ",props2," \"",sdfFile[i],"\"  >  \"",outputFile[i],"\"",sep=""),wait=TRUE)
}



