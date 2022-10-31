global.gates <- list()
global.gates[["Cancer"]] <- list(
              Pos=c('CK7','ECadherin'),
              Neg=c())
global.gates[["Lymphoid1"]] <- list(
              Pos=c('CD3d', 'CD8a'),
              Neg=c('CK7','ECadherin'))
global.gates[["Lymphoid2"]] <- list(
              Pos=c('CD3d', 'CD4'),
              Neg=c('CK7','ECadherin'))
global.gates[["Lymphoid3"]] <- list(
              Pos=c('CD20'),
              Neg=c('CK7','ECadherin'))
global.gates[["Myeloid1"]] <- list(
              Pos=c('CD11b','CD163'),
              Neg=c('CK7','ECadherin'))
global.gates[["Myeloid2"]] <- list(
              Pos=c('CD11b','CD15','CD68'),
              Neg=c('CK7','ECadherin'))
global.gates[["Myeloid3"]] <- list(
              Pos=c('CD11b','CD11c'),
              Neg=c('CK7','ECadherin'))
global.gates[["Stromal1"]] <- list(
              Pos=c('Vimentin','aSMA','CD31'),
              Neg=c('CK7','ECadherin'))
global.gates[["Stromal2"]] <- list(
              Pos=c('Vimentin','aSMA'),
              Neg=c('CK7','ECadherin'))
#global.gates[["Others"]] <- list(
#              Pos=c(),
#              Neg=unique(unlist(lapply(global.gates, function(x) x$Pos))))


immune.gates <- list()
immune.gates[["CD8.T.cells"]] <- list(
  Pos=c('CD8a','CD3d'),
  Neg=c('CD68','CD11c','CD163','FOXP3','CD4','CD15'))
immune.gates[["Lymphoids1"]] <- list(
  Pos=c('CD4','CD3d','FOXP3'),
  Neg=c('CD68','CD11c','CD163','CD15'))
immune.gates[["Lymphoids2"]] <- list(
  Pos=c('CD4','CD3d'),
  Neg=c('CD68','CD11c','CD163','CD15'))
immune.gates[["Myeloids1"]] <- list(
  Pos=c('CD11c','CD163'),
  Neg=c('CD3d'))
immune.gates[["Myeloids2"]] <- list(
  Pos=c('CD15','CD68'),
  Neg=c('CD3d'))
immune.gates[["Myeloids3"]] <- list(
  Pos=c('CD68','CD11c','CD163','CD15'),
  Neg=c('CD3d'))
immune.gates[["B cells"]] <- list(
  Pos=c('CD20'),
  Neg=c('CD68','CD11c','CD163','CD15','CD3d'))
immune.gates[["Other.immune"]] <- list(
  Pos=unique(unlist(lapply(immune.gates, function(x) x$Pos)[!grepl("Other.immune",names(immune.gates))])),
  Neg=c())



Lymphoid.gate <- list()
Lymphoid.gate[["CD4.T.cells"]] <- list(
  Pos=c('CD4'),
  Neg=c('FOXP3'))
Lymphoid.gate[["T.regs"]] <- list(
  Pos=c('FOXP3','CD4'),
  Neg=c())



Myeloid.gate <- list()
Myeloid.gate[["CD163.MP"]] <- list(
  Pos=c('CD163'),
  Neg=c('CD68','CD11c'))
Myeloid.gate[["CD68.MP"]] <- list(
  Pos=c('CD68'),
  Neg=c('CD11c','CD163'))
Myeloid.gate[["CD11c.MY"]] <- list(
  Pos=c('CD11c'),
  Neg=c('CD68','CD163'))
Myeloid.gate[["CD15.MY"]] <- list(
  Pos=c('CD15'),
  Neg=c('CD11c','CD163','CD68')) #CD68 was highly correlated with CD15
Myeloid.gate[["Other.MY"]] <- list(
   Pos=unique(unlist(lapply(Myeloid.gate, function(x) x$Pos)[!grepl("Other.MY",names(Myeloid.gate))])),
   Neg=c())


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



qc.immune.gates <- list()
qc.immune.gates[["CD8.T.cells"]] <- list(
  Pos=c('CD8a','CD3d'),
  Neg=c('CD11b','CD163','CD68','CD20'))
qc.immune.gates[["Lymphoids"]] <- list(
  Pos=c('CD4','CD3d','CD57'),
  Neg=c('CD11b','CD163','CD68','CD20','CD15'))
qc.immune.gates[["Myeloids1"]] <- list(
  Pos=c('CD11c','CD11b','CD207'),
  Neg=c('CD3d','CD20','CD4'))
qc.immune.gates[["Myeloids2"]] <- list(
  Pos=c('CD163','CD68','CD11b','CD15'),
  Neg=c('CD3d','CD20','CD4'))
qc.immune.gates[["B cells"]] <- list(
  Pos=c('CD20'),
  Neg=c('CD11b','CD163','CD68','CD3d'))
qc.immune.gates[["CancerImmune"]] <- list(
    Pos=c('CK7','ECadherin'),
    Neg=unique(c('CD57', 'CD15', unlist(lapply(immune.gates, function(x) x$Pos)[!grepl("CancerImmune",names(immune.gates)) & !grepl("CD8.T.cells",names(immune.gates)) & !grepl("Other.immune",names(immune.gates)) & !grepl("Lymphoids",names(immune.gates))]))))
qc.immune.gates[["Other.immune"]] <- list(
  Pos=unique(unlist(lapply(immune.gates, function(x) x$Pos)[!grepl("CancerImmune",names(immune.gates)) & !grepl("Other.immune",names(immune.gates)) & !grepl("Background_others",names(immune.gates))])),
  Neg=c())
