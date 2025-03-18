global.gates <- list()
# global.gates[["Cancer"]] <- list(
#               Pos=c('CK7','ECadherin'),
#               Neg=c())
global.gates[["Immune1"]] <- list(
              Pos=c('CD3d', 'CD8a'),
              Neg=c())
global.gates[["Immune2"]] <- list(
              Pos=c('CD3d', 'CD4'),
              Neg=c())
global.gates[["Immune3"]] <- list(
              Pos=c('CD3d', 'FOXP3'),
              Neg=c())
# global.gates[["B.cells"]] <- list(
#               Pos=c('CD20'),
#               Neg=c('aSMA'))
global.gates[["Immune4"]] <- list(
              Pos=c('CD68'),
              Neg=c())
global.gates[["Immune5"]] <- list(
              Pos=c('CD163'),
              Neg=c())
global.gates[["Immune6"]] <- list(
              Pos=c('CD207'),
              Neg=c())
global.gates[["Immune7"]] <- list(
              Pos=c('CD11c'),
              Neg=c())
global.gates[["Immune8"]] <- list(
              Pos=c('MHCII'),
              Neg=c('aSMA'))
global.gates[["Stromal1"]] <- list(
              Pos=c('Vimentin','aSMA'),
              Neg=c('CD3d','CD11c','CD163','CD207'))
global.gates[["Stromal2"]] <- list(
              Pos=c('Eccentricity'),
              Neg=c('CD3d','CD11c','CD163','CD207'))


immune.gates <- list()
immune.gates[["CD8.T.cells"]] <- list(
  Pos=c('CD8a','CD3d'),
  Neg=c('CD11c','CD163','FOXP3','CD4'))
immune.gates[["Lymphoids1"]] <- list(
  Pos=c('CD4','CD3d','FOXP3'),
  Neg=c('CD11c','CD163'))
immune.gates[["Lymphoids2"]] <- list(
  Pos=c('CD4','CD3d'),
  Neg=c('CD11c','CD163'))
immune.gates[["Myeloids1"]] <- list(
  Pos=c('CD163','CD68'),
  Neg=c('CD3d'))
immune.gates[["Myeloids2"]] <- list(
   Pos=c('CD11c'),
   Neg=c('CD3d'))
immune.gates[["Myeloids3"]] <- list(
  Pos=c('CD207'),
  Neg=c('CD3d'))
immune.gates[["Other.immune"]] <- list(
  Pos=unique(unlist(lapply(immune.gates, function(x) x$Pos)[!grepl("Other.immune",names(immune.gates))])),
  Neg=c())

Lymphoid.gate <- list()
Lymphoid.gate[["CD4.T.cells"]] <- list(
  Pos=c('CD4'),
  Neg=c('FOXP3','PD1'))
Lymphoid.gate[["T.regs"]] <- list(
  Pos=c('FOXP3','CD4'),
  Neg=c())
Lymphoid.gate[["CD4.CD45RO.T.cells"]] <- list(
  Pos=c('CD45RO'),
  Neg=c('FOXP3'))
Lymphoid.gate[["CD4.PD1.T.cells"]] <- list(
  Pos=c('PD1'),
  Neg=c('FOXP3'))

cd8.gates <- list()
cd8.gates[["CD8.PD1.T.cells"]] <- list(
  Pos=c('PD1'),
  Neg=c())
cd8.gates[["CD8.effector.T.cells"]] <- list(
  Pos=c('CD8a'),
  Neg=c('PD1'))
cd8.gates[["CD8.CD45RO.T.cells"]] <- list(
  Pos=c('CD45RO'),
  Neg=c())


Myeloid.gate <- list()
Myeloid.gate[["CD163.MP"]] <- list(
  Pos=c('CD163'),
  Neg=c('CD11c'))
Myeloid.gate[["CD207.MY"]] <- list(
  Pos=c('CD207'),
  Neg=c('CD163','CD11c'))
Myeloid.gate[["CD11c.MY"]] <- list(
  Pos=c('CD11c'),
  Neg=c('CD163'))
Myeloid.gate[["CD68.MP"]] <- list(
   Pos=c('CD68'),
   Neg=c('CD163','CD11c','CD207'))
Myeloid.gate[["Other.MY"]] <- list(
   Pos=c(),
   Neg=unique(unlist(lapply(Myeloid.gate, function(x) x$Pos)[!grepl("Other.MY",names(Myeloid.gate))])))
