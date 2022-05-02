##The following function combines the data from several columns of a 
##data frame into a single column. It then creates new columns that specify
##how the original columns were different from one another (i.e. factor columns).
##The output is a data frame with number of rows equal to the number of rows in 
##the original data frame times the number of columns combined together. The
##columns that were not combined together from the original data frame are left
##unchanged (i.e. they have the same column names and values - except that they
##will be repeated 'n' times, where 'n' is the number of columns combined).

##For example, if you had a data frame that had FPKM values in the columns, and
##wanted to combine the FPKM values into a single column and create factor
##columns specifying how the samples differ, this function would be used. 

##DataFrame should be the original data frame that has some columns you wish
##to combine together.

##MetaData is a data frame that specifies how each of the original columns
##differ from one another. There must be one column called "Column.Names" that 
##specifies the column names of DataFrame that you wish to combine together 
##into a single column. The remaining columns should be a set of values that 
##are unique to each original column name. The names of these columns of 
##MetaData will be set as column names in the new data frame that is created. 
##For example, if I had FPKM values of 2 replicates of 2 strains, my MetaData
##data frame should have two columns called "Replicate" and "Strain", which 
##specify the replicate number and strain name of each column name in the 
##Column.Names column of MetaData. Then, the new data frame will have two
##columns called "Replicate" and "Strain" that specify those values.

##ColNameCol should be a number specifying which column of the MetaData data 
##frame has the column names of the DataFrame data frame.

##NewColName will be set as the column name of the column in the new data frame
##containing all the columns combined into one from DataFrame. 

ColumnComb<-function(DataFrame,MetaData,ColNameCol,NewColName) {
  
  #Set the columns of DataFrame to be combined together and those to be left
  #alone. Also, set the columns of MetaData to become new columns in the new
  #data frame. 
  ColsToCombine<-match(MetaData[,ColNameCol],colnames(DataFrame))
  ColsToLeave<-(1:ncol(DataFrame))[-ColsToCombine]
  NewColsToAdd<-(1:ncol(MetaData))[-ColNameCol]
 
  #Collapse the columns to be combined into a single vector.
  UnlistData<-as.numeric(unlist(DataFrame[,ColsToCombine])) 
  
  #Create a new data frame whose columns correspond to the unchanged columns
  #of DataFrame, the columns to be combined from DataFrame and the factor 
  #columns specifying how those original columns combined were different (i.e. 
  #strain or replicate for FPKM values). 
  NewDataFrame<-as.data.frame(matrix(NA,nrow=length(UnlistData),ncol=ncol(DataFrame)-nrow(MetaData)+ncol(MetaData)),stringsAsFactors=FALSE)
  colnames(NewDataFrame)<-c(colnames(DataFrame)[ColsToLeave],NewColName,colnames(MetaData)[NewColsToAdd])
  
  #Add in the unchanged columns from DataFrame first (they need to be repeated
  #nrow(MetaData) times so that they are of appropriate length).
  for (k in ColsToLeave) {
    NewDataFrame[,colnames(NewDataFrame)==colnames(DataFrame)[k]]<-rep(DataFrame[,k],nrow(MetaData))
  }
  
  #Then, add in the new factor columns which specify how the combined columns
  #differ from one another (i.e. strain, replicate, condition, etc.)
  for (k in NewColsToAdd) {
    NewDataFrame[,colnames(NewDataFrame)==colnames(MetaData)[k]]<-rep(MetaData[,k],each=nrow(DataFrame))
  }
  
  #Add in the combined column (i.e. UnlistData)
  NewDataFrame[,colnames(NewDataFrame)==NewColName]<-UnlistData
  
  #Return the newly created NewDataFrame
  return(NewDataFrame)
}

