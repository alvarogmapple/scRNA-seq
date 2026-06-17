



library(CellChat)
library(Seurat)

library(shiny)
library(bsicons)
library(presto)

###1.Data input and preprocessing

nadim_seurat<-readRDS("C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/5.SINGLER (anotacion)/nadim_seurat_annotated.rds")

# Crear objeto cellchat

nadim_cellchat <- createCellChat(object = nadim_seurat, group.by = "SingleR_cluster_fine", assay = "RNA")


#Seleccionamos la base de datos a utilizar
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB) 
dplyr::glimpse(CellChatDB$interaction)  #te muestra todos tus datos


##Select all CellChatDB database
CellChatDB.use <- CellChatDB

#Añade como slot la DB a utilizar
nadim_cellchat@DB <- CellChatDB.use 

#Subsetea los datos de genes de señalización para ahorrar cpu
nadim_cellchat <- subsetData(nadim_cellchat)


#future hace computacion en paralelo, se divide el trabajo
#multisession abre varias sesiones, en las que 4 nucleos (workers) van a trabajar a la vez
#Esto amplia la memoria RAM utilizable, separandolo en varios planos
future::plan("multisession", workers = 4) 

#Para esta función hace falta paquete presto (test de wilcoxon rápido y optimizado)

#El rank rum test de wilcoxon compara la expresion de cada gen en un cluster con la de los otros clusters

#Overexpressed ligands or receptors (genes en general)
nadim_cellchat <- identifyOverExpressedGenes(nadim_cellchat)

#Overexpressed interaction (either ligand or receptor is overexpressed)
nadim_cellchat <- identifyOverExpressedInteractions(nadim_cellchat)




###2.Inference of cell–cell communication networks


#Inferimos ya la cell-cell comunication (esta es la default)
nadim_cellchat <- computeCommunProb(nadim_cellchat, type = "triMean", trim = NULL, 
                              raw.use = TRUE) 
    saveRDS(nadim_cellchat, file = "C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/6. CELLCHAT/1.nadim_cellchat_computeCommunProb.rds")
    #raw.use TRUE utiliza los datos raw, y FALSE los smooth

    nadim_cellchat<-readRDS("C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/6. CELLCHAT/1.nadim_cellchat_computeCommunProb.rds")

#Se   FILTRAN las comunicaciones, por DEAFAULT el nº de celulas minimo por grupo para que se considere comunicacion = 10, pero cambiar segun el numero de celulas posibles para comunicacion que tengamos disponibles
nadim_cellchat <- filterCommunication(nadim_cellchat, min.cells = 10)

#Cell-cell communication a nivel de pathway
nadim_cellchat <- computeCommunProbPathway(nadim_cellchat)


#(A) Perform calculation across all the cell groups
nadim_cellchat <- aggregateNet(nadim_cellchat)


#Exportar objeto cellchat
saveRDS(nadim_cellchat, file = "C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/6. CELLCHAT/2.nadim_cellchat_network.rds")

nadim_cellchat<-readRDS("C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/6. CELLCHAT/2.nadim_cellchat_network.rds")

###3.Visualization of cell–cell communication networks

#subsetCommunication: Crea data frame con todas las comunicaciones celulares inferidas
##! Esto da prob, luego count y weight lo dara aggregateNet, que con eso se general las figuras
COMMUNICATIONS<-subsetCommunication(nadim_cellchat)

#all the signaling pathways
 pathways.show.all <- nadim_cellchat@netP$pathways

#Select one pathway:
 pathways.show <- c("APP", "SPP1", "MIF", "FN1","COLLAGEN")


##(A) Circle plot
netVisual_aggregate(nadim_cellchat, signaling = pathways.show.all, layout 
                    = "circle", signaling.name="")

netVisual_aggregate(nadim_cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(nadim_cellchat, signaling = "APP", layout = "circle", signaling.name="APP")
netVisual_aggregate(nadim_cellchat, signaling = "SPP1", layout = "circle", signaling.name="SPP1")
netVisual_aggregate(nadim_cellchat, signaling = "MIF", layout = "circle", signaling.name="MIF")
netVisual_aggregate(nadim_cellchat, signaling = "FN1", layout = "circle", signaling.name="FN1")
netVisual_aggregate(nadim_cellchat, signaling = "COLLAGEN", layout = "circle", signaling.name="COLLAGEN") #Todas adipocytes (source)


netVisual_individual(nadim_cellchat, signaling = "APP",
                     pairLR.use = "APP_CD74", layout = "circle")

##(B) Hierarchy plot

vertex.receiver = seq(1,4) 
netVisual_aggregate(nadim_cellchat, signaling = pathways.show.all, signaling.name="", 
                    layout = "hierarchy", vertex.receiver = vertex.receiver, title= 2
)
#(C) Chord diagram

par(mfrow=c(1,1)) 
netVisual_aggregate(nadim_cellchat, signaling = pathways.show.all,  
                    layout = "chord", signaling.name = "")

netVisual_aggregate(nadim_cellchat, signaling = "SPP1", layout = "circle")

#(D) Heat map plot
par(mfrow=c(1,1)) 
netVisual_heatmap(nadim_cellchat, signaling = pathways.show.all,  
                  color.heatmap = "Reds")


#APLICACION INTERESANTE



#Todas las interacciones recibidas por un grupo celular:
###sources.use los SENDERS y target.use los receivers
netVisual_chord_gene(nadim_cellchat, sources.use = c(2,3,4), 
                     +                      targets.use = 1, legend.pos.x = 15)

#también se pueden comparar cell groups con una bubble plot:
netVisual_bubble(nadim_cellchat, sources.use = 4, targets.use = 
                   c(5:11), remove.isolate = FALSE)
####
### TOP 5 RUTAS INTERACTIONS bubble plot
   netVisual_bubble(nadim_cellchat, signaling = pathways.show, remove.isolate = FALSE)
######



###4.Systematic analysis of cell–cell communication

#Compute centrality scores of inferred-cell-cell communication network

  #WEIGHT: probabilidad de comunicación acumulada de todas esas interacciones
  #COUNT: número de interacciones ligando-receptor significativas que están activas en cada pathway
   
#RANKING pathways
   rankNet(nadim_cellchat, slot.name = "netP", mode = "single", 
        measure = "weight", sources.use = NULL, targets.use = NULL, 
        stacked = T, do.stat = FALSE)
   rankNet(nadim_cellchat, slot.name = "netP", mode = "single", 
           measure = "count", sources.use = NULL, targets.use = NULL, 
           stacked = T, do.stat = FALSE)


nadim_cellchat <- netAnalysis_computeCentrality(nadim_cellchat, slot.name = "netP")

      saveRDS(nadim_cellchat, file="C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/6. CELLCHAT/nadim_cellchat_centrality.rds")
      nadim_cellchat<-readRDS("C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/6. CELLCHAT/nadim_cellchat_centrality.rds")
      

 pathways.show<-c("COLLAGEN","APP","SPP1","FN1","LAMININ","MIF","PECAM1","MK","THBS")
 netAnalysis_signalingRole_network(nadim_cellchat, signaling = pathways.show)


#(B) Visualize dominant senders (sources) and receivers (targets) in a 2D space
 netAnalysis_signalingRole_scatter(nadim_cellchat)

#(C) Identify the major contributing signaling events of each cell group
#(i) Identify the major outgoing signaling events
 HM1 <- netAnalysis_signalingRole_heatmap(nadim_cellchat, pattern = 
                                           "outgoing",width = 10,
                                         height = 13) 

#(ii) Identify the major incoming signaling events
 HM2 <- netAnalysis_signalingRole_heatmap(nadim_cellchat, pattern = 
                                           "incoming",width = 10,
                                         height = 13
                                         ) 

#(iii) Show the major outgoing and incoming signaling events together
 HM1 + HM2
 #squares verdes son la relative signaling strength de la ruta en ese cell group
 #La barplot de arriba es la signaling strength total del cell group sumando todas las rutas
 #La barplot de la derecha es la signaling stregth total de la ruta sumando todos los cell groups (al reves de la de arriba)



#se agrupan las vias de señalización según su similaridad 
#en cuanto a RED DE COMUNICACION CELULAR (estructural o funcional)  

 library(reticulate)
 reticulate::py_install(packages = "umap-learn")
 
  ##(A) Functional similarity analysis
       #Patrones de quién emite y recibe¿Qué pathways tienen el mismo rol funcional en el TME?

#(i) Compute the functional similarity between any pair of inferred networks
nadim_cellchat <- computeNetSimilarity(nadim_cellchat, type = "functional")

#(ii) Perform manifold learning of inferred communication networks
nadim_cellchat <- netEmbedding(nadim_cellchat, type = "functional")

#(iii) Perform clustering of inferred communication networks
nadim_cellchat <- netClustering(nadim_cellchat, type = "functional")

#(iv) Visualize inferred communication networks in a 2D space
netVisual_embedding(nadim_cellchat, type = "functional",  
                    label.size = 3.5)
#(v) (Optional) Zoom in each group of signaling pathways in a 2D space
netVisual_embeddingZoomIn(nadim_cellchat, type = "functional",  
                          nCol = 2)


##(B) Structure similarity analysis
     #Topología de la red (qué nodos conecta)¿Qué pathways conectan los mismos tipos celulares?
nadim_cellchat <- computeNetSimilarity(nadim_cellchat, type = "structural") 

nadim_cellchat <- netEmbedding(nadim_cellchat, type = "structural") 

nadim_cellchat <- netClustering(nadim_cellchat, type = "structural") 

netVisual_embedding(nadim_cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(nadim_cellchat, type = "structural", nCol = 2)



computeNetSimilarityPairwise(nadim_cellchat, slot.name = "netP", c(1,2))

#rankSimilarity(nadim_cellchat)

saveRDS(nadim_cellchat, file="C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/6. CELLCHAT/4.nadim_cellchat_similarity.rds")

nadim_cellchat<-readRDS("C:/Users/alvar/OneDrive/Escritorio/BQ/MASTERº/TFM/6. CELLCHAT/4.nadim_cellchat_similarity.rds")

