#Importing the libs
library(ggplot2)
library(reshape2)
library(data.table)
library(cluster)
library(caTools) 
library(tidyverse)
library(factoextra)

#Setting the sys-path
root.dir = dirname(rstudioapi::getSourceEditorContext()$path)

#Import the data
data = readRDS(paste0(root.dir, "/DepMap19Q1_allData.RDS")) #read in the data (with paste0 you can concatinate strings; here some paths this makes your code very flexible)
length(data) #get some info about the data (the dimensions here)
lapply(data, dim) #get some impression about the data
lapply(data, function(a) head(a)) #look at the data 

#Extract the mutation data
mut <- data$mutation #pick the mutation data (this has a specific format and neets to be treated seperatelly)
'%!in%' <- function(x,y)!('%in%'(x,y)) #define an operator that will only pick the data that is NOT defined in the list; so the data that needs to be excluded so in our case: DNase, Methylation and RNA data should be excluded; so it will extract notin
dt_new <- lapply(which(names(data) %!in% "mutation"), function(a) data[[a]]) #extract only non-mutation data
names(dt_new) <- names(data)[which(names(data) %!in% "mutation")] #rename the data with the original names
length(dt_new) #look if it worked (now dim is 5 and NOT 6 (because mutation was removed))

#Define the samples you want to have a look at
sample_case = c("Breast Cancer") #define the samples you want to compare

#Extraktion der Mut Daten
samples = data$annotation$DepMap_ID[which(data$annotation$Primary.Disease == sample_case)]
ids = which(names(mut) %in% samples)
dat = lapply(ids, function(a) mut[[a]])

#Replace the colnames with real names
processed_data <- lapply(1:(length(dt_new)-1), function(a) {
  dat_picker <- dt_new[[a]] #pick one file at each iteration 
  col_names <- colnames(dat_picker) #extract the original column names
  annotation_new <- data$annotation[which(col_names %in% data$annotation$DepMap_ID),4] #get nice sample names in the order of the col_names list
  colnames(dat_picker) = annotation_new #rename the columns
  return(dat_picker)
})
names(processed_data) <- names(data)[which(names(data) %!in% c("mutation", "annotation"))] #rename the data
lapply(processed_data, dim) #look if the renaming worked
lapply(processed_data, head) #look if at the data

#Extract only kd.ceres values
BCCL_kd.ceres = processed_data$kd.ceres[,which(colnames(processed_data$kd.ceres) == sample_case)]
BCCL_Expression = processed_data$expression[,which(colnames(processed_data$expression) == sample_case)]
colnames(BCCL_Expression) <- samples
BCCL_Mutation = lapply(ids, function(a) mut[[a]])
names(BCCL_Mutation) = samples
BCCL_Annotation <- subset(data$annotation, Primary.Disease == sample_case)
tPatients_ID <- t(BCCL_Annotation[1]) 

#Initialize the statements
pre_mutImpact <- BCCL_kd.ceres
colnames(pre_mutImpact)<-tPatients_ID

#Initialize the empty matrix
mutImpact <- matrix(nrow = nrow(pre_mutImpact), ncol = ncol(pre_mutImpact))
colnames(mutImpact)<-colnames(pre_mutImpact)
rownames(mutImpact)<-rownames(pre_mutImpact)

#Add the data to the matrix
for (j in 1:ncol(pre_mutImpact)){
  lineID <- colnames(pre_mutImpact)[j] #select a column-name aka cell line
  for (i in 1:nrow(pre_mutImpact)){
    GOI <- rownames(pre_mutImpact)[i] #select a gene
    if (GOI %in% BCCL_Mutation[[lineID]]$Hugo_Symbol){
      mutImpact[i, j] <- pre_mutImpact[i, j]
    } else {
      mutImpact[i, j] <- NA
    }#replace CERES value with NA in case the gene is not mutated
  }
}

#Generate an input matrix for kMeans clustering
mutImpact_c = mutImpact[rowSums(is.na(mutImpact)) != ncol(mutImpact),]# remove all genes with all NA values
rownamesmutImpc = rownames(mutImpact_c)
mutImpact_c = cbind(rownamesmutImpc, mutImpact_c[,2:ncol(mutImpact_c)]) # insert new column with row names
rownameskdc2 = rownames(BCCL_kd.ceres)

BCCL_kd.ceres_3 = cbind(rownameskdc2, BCCL_kd.ceres[,2:ncol(BCCL_kd.ceres)])[cbind(rownameskdc2, BCCL_kd.ceres[,2:ncol(BCCL_kd.ceres)])[,"rownameskdc2"] %in%  mutImpact_c[,"rownamesmutImpc"],][,-c(1)]
BCCL_kd.ceres_3 = cbind(BCCL_kd.ceres_3, rowMeans(BCCL_kd.ceres_3))[rowMeans(cbind(BCCL_kd.ceres_3, rowMeans(BCCL_kd.ceres_3)))<= 0,] # deletion of genes with CERES > 0

km_CERES= kmeans(x = BCCL_kd.ceres_3, centers = 4, nstart = 10)
data_cluster = km_CERES$cluster
data_clus_3 = BCCL_kd.ceres_3[data_cluster == 3,] # without driver mutation label 

drivermut =  c("ERBB2","MYCBP", "PARP10", "PIK3CA") # create vector with driver mutation names

data_clus_3e = data_clus_3
data_clus_3e["Driver mutation"] <- sapply(1:nrow(data_clus_3), function(j){
  if(rownames(data_clus_3e)[j] %in% drivermut){
    data_clus_3e[j,1] = "driver mutation"
  } else {
    data_clus_3e[j,1] = "no driver mutation"
  }
}) # labeling driver mutations with "driver mutation" in data_clus_3e 

#Make a statistical test
mutImpact_r <- mutImpact[rowSums(is.na(mutImpact)) != ncol(mutImpact), ] #_r=reduced: get rid of all the ONLY NA rows (this will save computation)
mutImpact_kd.ceres <- BCCL_kd.ceres[rownames(BCCL_kd.ceres) %in% rownames(mutImpact_r),] # create matrix containing CERES scores of genes which are mutated once at minimum
colnames(mutImpact_kd.ceres) <- tPatients_ID 
testData <- lapply(seq_along(drivermut), function(a) {
  driverMutPicker <- drivermut[a] #pick a driver mutation
  driverMutData <- mutImpact_kd.ceres[driverMutPicker,] #get ceres scores of driver mutation
  
  outputData <- sapply(1:nrow(mutImpact_kd.ceres), function(b) {
    refGeneData <- mutImpact_kd.ceres[b,] #get ceres scores of reference gene (potential SST)
    
    if (rownames(refGeneData) != driverMutPicker) {
      out <- wilcox.test(as.numeric(driverMutData), as.numeric(refGeneData), paired = TRUE)$p.value #get p-value of wilcoxon signed-rank test
      out <- as.data.frame(out, rownames(mutImpact_r)[b]) 
      return(out)
    }
  })
  test = do.call(rbind, outputData)
  rownames(test) <- rownames(mutImpact_r)[which(rownames(mutImpact_r) != driverMutPicker)]
  return(test)
})
names(testData) <- drivermut
lapply(testData, function(a) head(a)) #look at the data

#Get the 2nd site targets
transformedData <- lapply(seq_along(testData), function(a) {
  data <- as.data.frame(testData[[a]], col.names = names(testData[[a]])) #create new data frame containing testData for ERBB2
  colnames(data) = drivermut[a] #rename column
  data <- data %>%rownames_to_column() %>% top_n(15, data[,1]) %>% column_to_rownames() #select 15 genes with highest p-value
  return(data)
})
names(transformedData) <- drivermut
lapply(transformedData, function(a) print(a))




####they still have to rewrite this a bit; I was running out of time:)
#Furhter process the data for a linear regression model
finalRegressionData <- lapply(seq_along(transformedData), function(a) {
  dtPicker <- transformedData[[a]]
  dtPicker$mutType <- sapply(seq_along(rownames(dtPicker)), function(b) {
    SSTPicker <- rownames(dtPicker)[b]
    mutTypeData <- c()
    for (i in 1:28) { #warum 28??? ncol????
      if (SSTPicker %in% BCCL_Mutation[[i]]$Hugo_Symbol) {
        out <- BCCL_Mutation[[i]]$Variant_Classification[a] 
        mutTypeData <- c(mutTypeData, out)
      }
    }
    dtPicker[b, 2] <- paste(mutTypeData, collapse = ", ")
  })
  return(dtPicker)
})
names(finalRegressionData) <- drivermut


#Make a linear regression model
mlr <- function(input_data, names, muts) {
  output <- lapply(seq_along(input_data), function(a){
    candidatas <- c(rownames(input_data[[a]]))
    dataPicker <- t(BCCL_kd.ceres[candidatas,])
    driverPicker <- t(BCCL_kd.ceres[muts[a], ])
    outData <- as.data.frame(cbind(dataPicker, driverPicker))
    return(outData)
  })
  names(output) <- names
  return(output)
}

test <- mlr(finalRegressionData, drivermut, driver_mut)
lapply(test, head)

lapply(seq_along(test), function(a) ggplot(test[[a]], aes(sample=test[[a]][,ncol(test[[a]])]))+stat_qq())

lapply(seq_along(input_data), function(a) {
  set.seed(123) #initialize the random numbers
  data <- input_data[[a]] #get the data
  colnames(data)[length(colnames(data))] <- "Predictor"
  
  split = sample.split(data[,ncol(data)], SplitRatio = 0.8) #split the dataset into 4/5 Training and 1/5 Testing dataset
  training_set = subset(data, split == TRUE) #use the labels to get the training data
  test_set = subset(data, split == FALSE) #dim(test_set) will give you know 10 --> 50/5*1 = 10; wuhu train/test split worked
  
  regressor = lm(formula = Predictor ~ ., #predict profit based on ALL (=.) the input variables for one company 
                 data = training_set)
  
  y_pred = predict(regressor, newdata = test_set) #predict the profit based on your testing data (this data the model did NEVER see and highly usefull to evaluat the performance)
  test_set$Prediction = y_pred #add your predictions to the dataset
  corVal <- cor.test(test_set$Predictor , test_set$Prediction, method = "spearman") #evaluate the model
  return(corVal)
})
