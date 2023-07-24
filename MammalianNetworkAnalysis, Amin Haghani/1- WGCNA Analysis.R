# Author: "Amin Haghani"
# Paper: "DNA Methylation Networks Underlying Mammalian Traits"
# Authors: "Amin Haghani1, 2 †*; Caesar Z. Li3, 4 †; Todd R. Robeck5; Joshua Zhang1; Ake T. Lu1, 2; Julia Ablaeva6; Victoria A. Acosta-Rodríguez7; Danielle M. Adams8; Abdulaziz N. Alagaili9, 10; Javier Almunia11; Ajoy Aloysius12; Nabil M.S. Amor13; Reza Ardehali14; Adriana Arneson15, 16; C. Scott Baker17; Gareth Banks18; Katherine Belov19; Nigel C. Bennett20; Peter Black21; Daniel T. Blumstein22, 23; Eleanor K. Bors17; Charles E. Breeze24; Robert T. Brooke25; Janine L. Brown26; Gerald Carter27; Alex Caulton28, 29; Julie M. Cavin30; Lisa Chakrabarti31; Ioulia Chatzistamou32; Andreas S. Chavez27, 33; Hao Chen34; Kaiyang Cheng35; Priscila Chiavellini36; Oi-Wa Choi37, 38; Shannon Clarke28; Joseph A. Cook39; Lisa N. Cooper40; Marie-Laurence Cossette41; Joanna Day42; Joseph DeYoung37, 38; Stacy Dirocco43; Christopher Dold44; Jonathan L. Dunnum39; Erin E. Ehmke45; Candice K. Emmons46; Stephan Emmrich6; Ebru Erbay47, 48, 49; Claire Erlacher-Reid43; Chris G. Faulkes50, 51; Zhe Fei3, 52; Steven H. Ferguson53, 54; Carrie J. Finno55; Jennifer E. Flower56; Jean-Michel Gaillard57; Eva Garde58; Livia Gerber59, 60; Vadim N. Gladyshev61; Rodolfo G. Goya36; Matthew J Grant62; Carla B. Green7; M. Bradley Hanson46; Daniel W. Hart20; Martin Haulena63; Kelsey Herrick64; Andrew N. Hogan65; Carolyn J. Hogg19; Timothy A. Hore66; Taosheng Huang67; Juan Carlos Izpisua Belmonte2; Anna J. Jasinska37, 68, 69; Gareth Jones70; Eve Jourdain71; Olga Kashpur72; Harold Katcher73; Etsuko Katsumata74; Vimala Kaza75; Hippokratis Kiaris76; Michael S. Kobor77; Pawel Kordowitzki78; William R. Koski79; Michael Krützen60; Soo Bin Kwon16, 15; Brenda Larison22, 80; Sang-Goo Lee61; Marianne Lehmann36; Jean-François Lemaître57; Andrew J. Levine81; Xinmin Li82; Cun Li83, 84; Andrea R. Lim1; David T.S. Lin85; Dana M. Lindemann43; Schuyler W.  Liphardt86; Thomas J. Little87; Nicholas Macoretta6; Dewey Maddox88; Craig O. Matkin89; Julie A. Mattison90; Matthew McClure91; June Mergl92; Jennifer J. Meudt93; Gisele A. Montano5; Khyobeni Mozhui94; Jason Munshi-South95; William J. Murphy96, 97; Asieh Naderi76; Martina Nagy98; Pritika Narayan62; Peter W. Nathanielsz83, 84; Ngoc B. Nguyen14; Christof Niehrs99, 100; Batsaikhan Nyamsuren101; Justine K. O'Brien42; Perrie O'Tierney Ginn72; Duncan T Odom102, 103; Alexander G. Ophir104; Steve Osborn105; Elaine A. Ostrander65; Kim M. Parsons46; Kimberly C. Paul81; Amy B. Pedersen87; Matteo Pellegrini106; Katharina J. Peters60, 107; Jessica L. Petersen108; Darren W. Pietersen109; Gabriela M. Pinho22; Jocelyn Plassais65; Jesse R. Poganik61; Natalia A. Prado110, 26; Pradeep Reddy111, 2; Benjamin Rey57; Beate R. Ritz112, 113, 81; Jooke Robbins114; Magdalena Rodriguez115; Jennifer Russell105; Elena Rydkina6; Lindsay L. Sailer104; Adam B. Salmon116; Akshay Sanghavi73; Kyle M. Schachtschneider117, 118, 119; Dennis Schmitt120; Todd Schmitt64; Lars Schomacher99; Lawrence B. Schook117, 121; Karen E. Sears22; Ashley W. Seifert12; Aaron B.A. Shafer122; Anastasia V. Shindyapina61; Melanie Simmons45; Kavita Singh123; Ishani Sinha22; Jesse Slone67; Russel G. Snell62; Elham Soltanmohammadi76; Matthew L. Spangler108; Maria Spriggs21; Lydia Staggs43; Nancy Stedman21; Karen J. Steinman124; Donald T Stewart125; Victoria J. Sugrue66; Balazs Szladovits126; Joseph S. Takahashi7, 127; Masaki Takasugi6; Emma C. Teeling128; Michael J. Thompson106; Bill Van Bonn129; Sonja C. Vernes130, 131; Diego Villar132; Harry V. Vinters133; Ha Vu15, 16; Mary C. Wallingford72; Nan Wang37, 38; Gerald S. Wilkinson8; Robert W. Williams134; Qi Yan3, 2; Mingjia Yao3; Brent G. Young54; Bohan Zhang61; Zhihui Zhang6; Yang Zhao6; Peng Zhao14, 135; Wanding Zhou136, 137; Joseph A. Zoller3; Jason Ernst15, 16; Andrei Seluanov138; Vera Gorbunova138; X. William Yang37, 38; Ken Raj139; Steve Horvath1, 2 *"
# Date: "07-24-2023"
# output: html_document

BiocManager::install("AnnotationDbi")
install.packages("WGCNA_1.69-81.tar.gz",repos=NULL, type="source")

library(WGCNA)

bVales <- readRDS("DNAmData.RDS")
datExpr <- t(bVales)
rm(bVales)

powers = c(c(1:10), seq(from = 10, to=160, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")


sizeGrWindow(12, 9)
pdf(file = "eigengene network merged.pdf", wi = 9, he = 6)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
dev.off()


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$fitIndices


softPower = 12
adjacency = adjacency(datExpr, power = softPower, type="signed");
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM
rm(adjacency)
rm(TOM)

geneTree = hclust(as.dist(dissTOM), method = "average")

minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)

# Test for module overlap and merge collinear modules (correlation of r>=0.85)
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.15
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

saveRDS(MEs, "Module MEs dynamic cut.RDS")
saveRDS(mergedMEs, "Module MEs merged.RDS")

dev.off()
sizeGrWindow(12, 9)
pdf(file = "geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic colors","Merged modules"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Dynamic colors
moduleColors = dynamicColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
colorh1 <- moduleColors
datME=moduleEigengenes(datExpr,colorh1, softPower = 12)$eigengenes
#Calculating KME
datKME=signedKME(datExpr, datME, outputColumnName="MM.")
datKME$modules <- colorh1

saveRDS(datKME, "dynamic cut module KME.RDS")

# Merged colors
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
colorh1 <- mergedColors
datME2=moduleEigengenes(datExpr,colorh1, softPower = 12)$eigengenes
#Calculating KME
datKME2=signedKME(datExpr, datME2, outputColumnName="MM.")
datKME2$modules <- colorh1

saveRDS(datKME2, "Merged modules KME.RDS")


MEs=orderMEs(datME)
sizeGrWindow(12, 9)
pdf(file = "eigengene network dynamic cut.pdf", wi = 9, he = 6)
plotEigengeneNetworks(MEs, "", marDendro = c(0,3.5,1,3.5), marHeatmap = c(3,4,1,2))
dev.off()


MEs=orderMEs(datME2)
sizeGrWindow(12, 9)
pdf(file = "eigengene network merged.pdf", wi = 9, he = 6)
plotEigengeneNetworks(MEs, "", marDendro = c(0,3.5,1,3.5), marHeatmap = c(3,4,1,2))
dev.off()





