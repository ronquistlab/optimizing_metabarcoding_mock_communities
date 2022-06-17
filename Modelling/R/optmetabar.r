# Help-function to read in results from Birch
# Optimal metabarcoding paper - June 2022
########################################################
# example use:
# test_H<-OptmetabarResults("test-H.json")

# LogZ: column 5
# k: column 6
# theta: column 7
# View for example logZ: boxplot(as.numeric(test_H[,5]))
########################################################

library(rjson)
OptmetabarResults<-function(filename) {
  data1<-fromJSON(file=filename)
  data2 <- do.call(rbind, data1)
  data3 <- do.call(rbind, data2[,6])
  final_data<-cbind(data2[,1:5],data3)
  return(final_data)
}
