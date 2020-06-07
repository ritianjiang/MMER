library(magrittr)
library(stringr)

oneStepTransMat<-function(st1,nd2){ #Generate 1-step transform matrix by two chr vectors
  return((table(st1,nd2) / rowSums(table(st1,nd2))) %>% as.matrix())
}

generateTransChain<-function(nucleotideSeqs){ #Generate chains for different start character
  split_Seq<-str_split(nucleotideSeqs,pattern="",simplify=T) #split sequences into matrix  
  length_Of_Seq<-ncol(split_Seq) #The length of each sequence. They should have same lengths.
  chain_File<-list()
  for(i in 1:(length_Of_Seq-1)){
    chain_File[[i]]<-oneStepTransMat(split_Seq[,i],split_Seq[,i+1])
  }
  (table(split_Seq[,1])/nrow(split_Seq))->freq_Init;freq_Init<-as.vector(freq_Init);
  names(freq_Init)<-names(table(split_Seq[,1])) #calculate the freqs of initial nucleotides
  chain_File<-c(list(freq_Init),chain_File)
  #names(chain_File)<-c("freqInit",paste0("TM",1:(length_Of_Seq-1)))
  return(chain_File)
}

predictByChain(chain_File,sequence){ #predict one sequence by chain-file
  split_Seq_One<-str_split(sequence,pattern="",simplify = T)
  if(length(split_Seq_One)!=(length(chain_File))){
    cat("The sequence's length doesn't match the chain_File")
    return()
  }
  mat_Elements<-sapply(2:length(split_Seq_One), function(x){ #get the Trans prob in each step
    return(chain_File[[x]][split_Seq_One[x-1],split_Seq_One[x]])
  }) 
  return(prod(mat_Elements) * chain_File[[1]][split_Seq_One[1]] %>% as.numeric())
}

plotChain<-function(chain_File,fromNuc,toNuc){ #visualize the prob chain of "one nucleotide to another"
  plot_Vec<-rep(0,length(chain_File)-1);
  for(i in 1:length(plot_Vec)){
    plot_Vec[i]<-chain_File[[i+1]][fromNuc,toNuc]
  }
  plot(1:length(plot_Vec),plot_Vec,"b",ylim=c(0,1))
}

#################Test#################


