
slopes_mat[SpeciesMat$MammalNumberHorvath=="10.1.3" & SpeciesMat$Tissue =="Skin", 11:13] = 
  NA


slopes_mat[SpeciesMat$MammalNumberHorvath=="8.7.9", 1] = 
  2*slopes_mat[SpeciesMat$MammalNumberHorvath=="8.7.9", 2]

slopes_mat[SpeciesMat$MammalNumberHorvath=="8.7.10",11] = NA
slopes_mat[SpeciesMat$MammalNumberHorvath=="8.7.4",11] = 
  slopes_mat[SpeciesMat$MammalNumberHorvath=="8.7.4",12] 

slopes_mat[SpeciesMat$MammalNumberHorvath=="8.9.1",2] = NA

slopes_mat[SpeciesMat$MammalNumberHorvath=="5.5.1",13] = 
  slopes_mat[SpeciesMat$MammalNumberHorvath=="5.5.1",11]

slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.1" & SpeciesMat$Tissue =="SVZ",12:13] = 
  slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.1" & SpeciesMat$Tissue =="SVZ",11]

slopes_mat[SpeciesMat$MammalNumberHorvath=="8.7.8", 1] = NA

if(len1 == "BivProm1+") {
  slopes_mat[SpeciesMat$MammalNumberHorvath=="9.5.10", 1] = NA
  slopes_mat[SpeciesMat$MammalNumberHorvath=="2.1.1", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="2.1.1", 2]*2
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="1.4.2", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="1.4.2", 2]*2
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="17.1.3", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="17.1.3", 2]*2
  
}
if(len1 == "BivProm2+") {
  slopes_mat[SpeciesMat$MammalNumberHorvath=="8.5.6", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="8.5.6", 2]*2
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="8.17.2", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="8.17.2", 2] = NA
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="5.2.1", 1:2] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="5.2.1", 3] *2
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.1" & SpeciesMat$Tissue =="Kidney",11] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.1" & SpeciesMat$Tissue =="Kidney",12]
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.21"&
               SpeciesMat$Tissue=="Liver", 12:13 ] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.21"&
                 SpeciesMat$Tissue=="Liver", 11 ]
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="21.2.1", 1:2] = NA
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="11.1.3", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="11.1.3", 2]
  
  
  slopes_mat[SpeciesMat$MammalNumberHorvath %in% c("1.3.3","1.3.9") , 11:13] = NA
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="1.2.3" & 
               SpeciesMat$Tissue=="Blood&Skin", 11:13] = NA
  slopes_mat[SpeciesMat$MammalNumberHorvath=="6.1.2", 11:13] = NA
  
  slopes_mat[SpeciesMat$MammalNumberHorvath %in% c("4.13.3","4.13.10"), 11:13] = NA
  
  #Lagenorhynchus obliquidens
  slopes_mat[SpeciesMat$MammalNumberHorvath=="4.13.4" & SpeciesMat$Tissue=="Blood", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="4.13.4" & SpeciesMat$Tissue=="Blood", 2]*2
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="1.1.1" & SpeciesMat$Tissue=="Skin", 1:3] = NA
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="1.1.1" & SpeciesMat$Tissue=="Blood", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="1.1.1" & SpeciesMat$Tissue=="Blood", 2]*2
  
}

if(len1 == "ReprPC1+") {
  slopes_mat[SpeciesMat$MammalNumberHorvath=="4.13.4" & SpeciesMat$Tissue=="Blood", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="4.13.4" & SpeciesMat$Tissue=="Blood", 2]*2 
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="9.1.3" & SpeciesMat$Tissue=="Liver", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="9.1.3" & SpeciesMat$Tissue=="Liver", 2]*2 
  
}


if(len1 == "BivProm3+") {
  slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.3" & SpeciesMat$Tissue =="Brain", 
             c(1:6, 11:13)] = NA
  
  slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.3" & SpeciesMat$Tissue =="Liver", 1] = 
    slopes_mat[SpeciesMat$MammalNumberHorvath=="9.9.3" & SpeciesMat$Tissue =="Liver", 2]
  
}