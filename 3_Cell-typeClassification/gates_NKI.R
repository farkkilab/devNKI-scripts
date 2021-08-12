global.gates <- list()
global.gates[["Tumor"]] <- list(
              Pos=c('CK7','ECadherin'),
              Neg=c())
global.gates[["Lymphoid1"]] <- list(
              Pos=c('CD3d', 'CD20', 'CD4'),
              Neg=c('CK7','ECadherin','Vimentin'))
global.gates[["Lymphoid2"]] <- list(
              Pos=c('CD3d', 'CD8a'),
              Neg=c('CK7','ECadherin','Vimentin'))
global.gates[["Myeloid cells"]] <- list(
              Pos=c('CD11c','MHCII','CD11b','CD207','CD163'),
              Neg=c('CK7','ECadherin', 'CD3d', 'CD20','Vimentin'))
global.gates[["Other"]] <- list(
            Pos=unique(unlist(lapply(global.gates, function(x) x$Pos))),
            Neg=c())
global.gates[["Stromal cells"]] <- list(
            Pos=c('Vimentin','aSMA'),
            Neg=c('CD3d', 'CD20', 'CD11c','CD163'))
#global.gates[["Background"]] <- list(
#            Pos=c('BG488', 'BG555', 'BG647'),
#            Neg=unique(unlist(lapply(global.gates, function(x) x$Pos))))
            #Neg=unique(unlist(lapply(global.gates, function(x) x$Pos))))

immune.gates <- list()
immune.gates[["CD8 T cells"]] <- list(
  Pos=c('CD8a','CD3d'),
  Neg=c('CD4','CD20','CD163'))
immune.gates[["CD4 T cells"]] <- list(
  Pos=c('CD4','CD3d'),
  Neg=c('CD8a','CD163','CD20'))
immune.gates[["Macrophages"]] <- list(
  Pos=c('IBA1','CD11b','CD163'),
  Neg=c('CD20','CD8a','CD3d','CD57','CD207','CK7','Ecadherin','PAX8'))
immune.gates[["B cells"]] <- list(
  Pos=c('CD20'),
  Neg=c('IBA1','CD4','CD8a','CD3d','CD163','CD1c','CD207','CK7','Ecadherin','PAX8'))
immune.gates[["Antigen presenting cells"]] <- list(
  Pos=c('CD207','CD1c'),
  Neg=c('IBA1','CD4','CD20','CD8a','CD57','CD3d','CD163','CK7','Ecadherin','PAX8'))
immune.gates[["NK cells"]] <- list(
  Pos=c('CD57'),
  Neg=c('PD1','IBA1','CD4','CD20','CD163','CD8a','CD1c','CD207','CK7','Ecadherin','PAX8'))
immune.gates[["Neutrophils"]] <- list(
  Pos=c('CD15'),
  Neg=c('IBA1','CD4','CD20','CD8a','CD57','CD3d','CD163','CD207','CK7','CD11b','CD1c','Ecadherin','PAX8'))
immune.gates[["Other"]] <- list(
  Pos=c(),
  Neg=c('IBA1','CD4','CD20','CD8a','CD57','CD163','CD11b','CD15','CD3d', 'CD1c'))

cd8.gates <- list()
cd8.gates[["Exhausted CD8 T cells"]] <- list(
  Pos=c('PD1','CD8a','CD3d'),
  Neg=c())
cd8.gates[["CD8 effector T cells"]] <- list(
  Pos=c('CD8a','CD3d'),
  Neg=c('PD1'))

cd4.gates <- list()
cd4.gates[["CD4 effector T cells"]] <- list(
  Pos=c('CD4','CD3d'),
  Neg=c('FOXP3'))
cd4.gates[["T regulatory cells"]] <- list(
  Pos=c('FOXP3','CD4','CD3d'),
  Neg=c())

macrophage.gates <- list()
macrophage.gates[["IBA1 CD11b Macrophages"]] <- list(
  Pos=c('CD11b','IBA1'),
  Neg=c('CD163'))
macrophage.gates[["CD163 Macrophages"]] <- list(
  Pos=c('CD163','IBA1'),
  Neg=c('CD11b'))
