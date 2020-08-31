# This code snippet accompanies the manuscript "Physicochemical and structural parameters contributing to the 
# antibacterial activity and efflux susceptibility of small molecule inhibitors of Escherichia coli", by 
# El Zahed et al. (2020).  It uses example data found at github.com/sfrench007/serf, and builds a random forest 
# model to predict molecules that are actively pumped in E. coli.  
#
# The code borrows from the elegant work by Richter et al. (2017; https://www.nature.com/articles/nature22308),
# which is found at github.com/HergenrotherLab/GramNegAccum.  

options(stringsAsFactors=FALSE)
library(caret)
library(doParallel)
library(ggplot2)
library(pROC)
library(GGally)

# Define the number of cores to use (typically ncores-1)
# (outfile="" is needed to return verbose information when tuning model in parallel)
cl <- makePSOCKcluster(7, outfile="") # ie. for 8 cores (4 physical, 4 virtual) ideal to use 7
registerDoParallel(cl)

# Reproducibility if required for publication purposes - this will save the seed to reproduce your data
use_seed <- sample(1:100000,1)  # Generate a random seed between 1-100000
set.seed(use_seed)

# Define a filename prefix for model and plots
prefix <- paste0("serf_seed",use_seed)  # Can be anything, defaulting here to include the seed used

# Load example datasets
training <- read.csv("examples/training_set.csv",sep="\t",header=TRUE) # N.B.: Near-zero variance already removed
test_set <- read.csv("examples/test_set.csv",sep="\t",header=TRUE)

# Subset the activity and properties
activity <- as.factor(training[,1])
properties <- training[,2:ncol(training)]
fullSet <- names(training)[names(training)!="Class"]

# Remove Near Zero variance and number of conformers if it's in there
isNZV <- nearZeroVar(training[,fullSet],saveMetrics=TRUE,freqCut=floor(nrow(training)/5))
fullSet <- rownames(subset(isNZV, !nzv))

# Some properties are extremely correlated, which complicates RF model building
cor_matrix <- cor(apply(training[,fullSet],1:2,as.numeric))
correlated <- findCorrelation(cor_matrix,0.95) # Use a 95% similarity
fullSet <- fullSet[-correlated]

# Define training control function (iterations, repeats, class summary, etc)
# A random search is used here over a grid-based search, with a high number of iterations
# and repeats.  A grid based search may be more robust, depending on the starting values
# used, but we prefer a randomized searching approach
ctrl <- trainControl(method= "repeatedcv",
    number=20,
    repeats=20,
    summaryFunction=twoClassSummary,
    classProbs=TRUE,
    savePredictions=TRUE,
    returnResamp="all",
    search="random",
    verboseIter=TRUE,
    allowParallel=TRUE)

# Number of properties possible at each tree node
mtryValues <- c(2,4,6,8,10,20) # This may vary, likely to be optimal at lower values

# Format the data and remove factors
finalSet <- data.frame(lapply(training[,fullSet],as.numeric))
colnames(finalSet) <- fullSet
plotThis <- data.frame(cbind(activity,data.frame(lapply(finalSet,as.numeric))))
colnames(plotThis)[2:length(plotThis)] <- fullSet

# Build the random forest model using caret with the rf method.  Pre-processing is centering and scaling
# but a PCA pre-processing step has also shown to be effective in past work.  
rf_model <- train(x=finalSet,
    y=factor(activity),
    method="rf",
    preProcess=c("center","scale"),
    ntree=2000,
    tuneGrid=data.frame(mtry=mtryValues),
    importance=TRUE,
    metric="ROC",
    trControl=ctrl)

# Calculate the confusion matrix, then save the model and confusion matrix
conf_matrix <- confusionMatrix(rf_model,norm = "none")
save(rf_model,file=paste(prefix,"_RF.Rdata",sep=""))
save(conf_matrix,file=paste(prefix,"_CM.Rdata",sep=""))

# Calculate important variables
nvar <- 10
all_props <- varImp(rf_model)
top_props <- all_props$importance[with(all_props$importance, order(-Active)),]
top_names <- rownames(top_props[1:nvar,]) 
final <- data.frame(name=rownames(top_props[1:nvar,]),importance=top_props$Active[1:nvar])
final <- data.matrix(top_props$Active[1:nvar])
rownames(final) <- rownames(top_props[1:nvar,])

rf_roc <- roc(response=rf_model$pred$obs,
    predictor=rf_model$pred$Active,
    levels=rev(levels(rf_model$pred$obs)))

# Plot tuning parameters
pdf(paste(prefix,"_tuning.pdf",sep=""),width=6,height=6)
plot(rf_model)
dev.off()

# Plot ROC
pdf(paste(prefix,"_ROC.pdf",sep=""),width=6,height=6)
plot(rf_roc,type="s",print.thres=c(.5),
     print.thres.pch=3,legacy.axes=TRUE,print.thres.pattern="",
     print.thres.cex=1.2,
     col="red",print.thres.col="red",
     print.auc=TRUE,print.auc.x=0.8,print.auc.y=0.6)
dev.off()

# Plot top <nvar> variables
pdf(paste(prefix,"_keyprops.pdf",sep=""),width=6,height=6)
dotchart(final)
dev.off()

# Plot matrix of important variable interactions
pdf(paste(prefix,"_matrix.pdf",sep=""),width = 14,height =14)
ggpairs(plotThis,mapping=aes(color=activity,alpha=0.5),columns=top_names,
    upper=list(continuous=wrap("cor",size=4,alignPercent=1)),
    lower=list(continuous=wrap("points",size=0.6)),
    diag=list(continuous="densityDiag"),axisLabels="show",
    title="Key Physicochemical Properties from Random Forest") +
    theme_linedraw(base_size=8) +
    theme(plot.title=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text=element_text(size=8),
        legend.position="top",
        legend.title=element_blank())
dev.off()

# Now predict from the test set, using the generated model
efflux_predictions <- predict(rf_model,newdata=test_set,type="prob")

