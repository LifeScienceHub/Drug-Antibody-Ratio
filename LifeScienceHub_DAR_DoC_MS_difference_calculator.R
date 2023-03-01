# ----- Packages:
library("readxl")
library("dplyr")
library("xlsx")
library("ggplot2")
library("ggpubr")

# load the MS data
# data containing MS spectra with proteoforms found and their intensity values - relative intensity values in % can be usually found on the spectra

df_MS<-read_excel("LifeScienceHub_example-data_intensities_deconvoluted_MS_analysis.xlsx", sheet= 2, col_names = TRUE, col_types = NULL, na = "", skip = 0)

####### Function for probable additional masses
# mz_range mass +/- value (value between 5-10 is reasonable for masses >5kDa)
# df_MS is addiotional, conjugated fragments to the protein/antibody e.g. drug or chelator


function_masses <- function(mz_range, df_MS){
  
  Bn_NOTA <- 581
  two_Bn_NOTA <- 1162
  three_Bn_NOTA <-1743
  Bn_NOTA_Cu <- 643
  two_Bn_NOTA_Cu<-1286
  three_Bn_NOTA_Cu<-1929
  
  fragments<-c(Bn_NOTA, two_Bn_NOTA, three_Bn_NOTA, Bn_NOTA_Cu, two_Bn_NOTA_Cu, three_Bn_NOTA_Cu)
  fragment_type<-c("Bn_NOTA", "two_Bn_NOTA", "three_Bn_NOTA", "Bn_NOTA_Cu", "two_Bn_NOTA_Cu", "three_Bn_NOTA_Cu")
  data_names<-c("mass", "fragment type")
  
  
  # mz_range limits (plus-minus m/z):
  # ------ Loading basic Data frames:
  masses<- df_MS$Mass
  intensities<-df_MS$Intensity
  
  
  ### creating evaluation table and adding value:
  data_names<-c("fragment type", "fragment", "mass1", "mass2", "delta", "mz_range")
  data<-data.frame(matrix(ncol=length(data_names), nrow=0))
  colnames(data)<-data_names
  #data[1, ]<-1:length(data_names) # first nmbr column, second row
  
  # ------------------- Listing of all Mass differences
  length_fragments<-length(fragments)
  length_fragmenttype<-length(data_names)
  
  data_fragments<-data.frame(matrix(ncol=length_fragmenttype, nrow=length_fragments))
  colnames(data_fragments)<-data_names
  data_fragments$mass<-fragments
  data_fragments$`fragment type`<-fragment_type
  print(data_fragments)
  
  
  ### -------------------------------- Function for data Analysis of all proteoforms with listing
  # function evaluates all proteoforms + the additional conjugated mass defined above
  k <- 1
  for (f in 1:length_fragments){
    frag<-data_fragments$mass[f]
    frag_type<-data_fragments$`fragment type`[f]
    f_upper_limit<-frag+mz_range
    f_lower_limit<-frag-mz_range
    for (i in 1:(length(masses)-1)){
      for(j in (i+1):length(masses)){
        v1 <- masses[i]
        v2 <- masses[j]
        difference<-abs(v1 - v2)
        diff1<-between(difference, f_lower_limit, f_upper_limit)
        if(diff1==TRUE){
          data[k, ]<-c(frag_type, frag, v1, v2, difference, mz_range)
          k<-k+1
        }
      }
    }
  }
  
  newdata<-data[order(data$mass1),]
  print(newdata)
  
  length_masses<-length(newdata$mass1)
  mass1<-newdata$mass1
  mass2<-newdata$mass2
  masses<- df_MS$Mass
  intensities<-df_MS$Intensity
  length_intensities<-length(intensities)
  
  for(i in 1:length_masses){
    
    for(k in 1:length_intensities){
      mass_compare_daughter_frag<-masses[k]
      mass_daughter_frag1<-mass1[i]
      mass_daughter_frag2<-mass2[i]
      #print(mass_compare_daughter_frag)
      #newdata$Added_Column<-ifelse(mass_compare_daughter_frag==mass_daughter_frag1, intensities[k], "NA")
      
      if (mass_compare_daughter_frag==mass_daughter_frag1){
        newdata$Mass1_intensity[i]<-intensities[k]
      }
      if (mass_compare_daughter_frag==mass_daughter_frag2){
        newdata$Mass2_intensity[i]<-intensities[k]
      }
    }
  }
  return(newdata)
}

data_evaluation <- function_masses(5, df_MS)


# print/export evaluation data file to elucidate the final DAR or DoC
# %-intensities should be taken from the spectra
# !!! careful, doupletts can be present as 2x Bn_NOTA and 1x Two_NOTA have the same MS difference to one parent proteoform, verify the data to exclude them for your final elucidation. 

print(data_evaluation)
write.csv(data_evaluation, "LifeScienceHub_example_MS_evaluation-file_range5.csv")





