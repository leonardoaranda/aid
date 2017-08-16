library(ggplot2)
library(car)
library(ggbiplot)
library(cluster)
library(MASS)
library(klaR)
library(AUC)
library(aplpack)
library(graphics)
library(pryr)
library(scatterplot3d)
library(caret)
library(forecast)
library(corrplot)
library(RColorBrewer)
library(gridExtra)

set.seed(83089)

par(ps = 8, cex = 1, cex.main = 1)

df <- read.csv("data/sites.csv",row.names=16)
df <- df[complete.cases(df),]
df$time <- as.numeric(df$time)
df <- subset(df,bounce > 5)
df <- subset(df,ppv < 20)
df$position <- as.factor(ifelse(df$global < 50000,1,0))


table(df$position,df$category)
cc <- c("red","blue")
ccc <- c("#6495ED")

# Preparo variables
vars <- c("bounce","paid_search","ppv","time","traffic_direct","traffic_display","traffic_search","traffic_mail","traffic_social","traffic_referrals")
df.num <- df[vars]
df.num <- as.data.frame(scale(df.num))
df.num.nos <- df[vars]

# Exploración gráfica

plots.pairs %<a-% pairs(df.num.nos,upper.panel = NULL,col=cc[as.numeric(df$position)],pch=20)
tables.means <- t(aggregate(df.num.nos,list(df$position),function(x) c(mean = round(mean(x, na.rm=TRUE), 2))))


# Análisis de tiempo en página

plots.time %<a-% hist(df$time,col="slateblue4",xlab="Tiempo de Visita",ylab="Frecuencia absoluta",main="Distribución del Tiempo de Visita", cex.main=0.8,cex.axis=0.8,cex.lab=0.8)
plots.timehist %<a-% boxplot(df$time ~ df$position,col="slateblue4", main="Tiempo de visita según posición",cex.main=0.8,cex.lab=0.8,cex.axis=0.8)
plots.time
plots.timehist

plots.qq1 %<a-% qqnorm(df$time, cex.main=0.8,cex.axis=0.8,cex.lab=0.8)
plots.qq2 %<a-% qqline(df$time, col = 2)
plots.qq1
plots.qq2

tables.shapiro <- shapiro.test(df$time)

tables.time <- aggregate(df$time,list(df$position),function(x){y <- shapiro.test(x); return(y$p.value)})
names(tables.time) <- c("position","p-value")

# Box-cox transformation

powerr <- powerTransform(time ~ position,data=df)
transs <- bcPower(df$time,powerr$lambda)

boxcox <-  aggregate(transs,list(df$position),function(x){y <- shapiro.test(x); return(y$p.value)})
plots.boxcox %<a-% boxplot(transs ~ df$position,col="slateblue4",main="Transformación Box-Cox",cex.main=0.8,cex.lab=0.8,cex.axis=0.8)
names(boxcox) <- c("position","p-value")



# Prueba de la mediana

tables.wilcox <- wilcox.test(time ~ position,data=df)

# Análisis del rebote

powerr <- powerTransform(time ~ position,data=df)
transs <- bcPower(df$time,powerr$lambda)

aggregate(transs,list(df$position),function(x){y <- shapiro.test(x); return(y$p.value)})


# Componentes principales

png(filename="fig/corr.png")
corrplot(cor(df.num),order="hclust",tl.col = "black",col=brewer.pal(10,"Spectral"))
dev.off()

df.num.pca <- prcomp(df.num)


plots.scree %<a-% screeplot(df.num.pca,type="lines",main="Scree Plot",cex.main=0.8,cex.axis=0.8,cex.lab=0.8)

comp1 <- df.num.pca$rotation[,1]
comp2 <- df.num.pca$rotation[,2]
comp3 <- df.num.pca$rotation[,3]

df$pca1 <- df.num.pca$x[,1]
df$pca2 <- df.num.pca$x[,2]
df$pca3 <- df.num.pca$x[,3]


plots.pc1 %<a-% barplot(comp1,las=2,main="Componente Principal #1", cex.names=0.7, cex.main=0.8,cex.axis=0.8,cex.lab=0.8,col="slateblue4") # Engagement
plots.pc2 %<a-% barplot(comp2,las=2,main="Componente Principal #2",cex.names=0.7, cex.main=0.8,cex.axis=0.8,cex.lab=0.8,col="slateblue4") # Awareness
plots.pc3 %<a-% barplot(comp3,las=2,main="Componente Principal #3",cex.names=0.7, cex.main=0.8,cex.axis=0.8,cex.lab=0.8,col="slateblue4") # Paid mkt


# Clustering


df.dist <- dist(df.num)
df.clust <- hclust(df.dist,method="average")
df.cof <- cophenetic(df.clust)
dist <- cor(df.dist,df.cof)

png(filename="fig/dend.png")
plot(as.dendrogram(df.clust),main="Clustering jerárquico",cex.main=0.8)
dev.off()


plot(as.dendrogram(hclust(df.dist,method="average")),main="Clustering jerárquico",cex.main=1)

df$cluster <- as.factor(kmeans(df.num,centers=4)$cluster)


plots.biplot <- ggbiplot(df.num.pca,groups=df$position,alpha=0.1) + theme_grey(base_size = 8) + geom_point(aes(colour=df$position),size=3,alpha=0.1) + scale_color_manual(values=cc) + labs(title="Position sobre componentes principales") + xlim(c(-5,5))  + ylim(c(-7,7)) + theme(panel.background = element_blank())

plots.biplot2 <- ggbiplot(df.num.pca,groups=df$cluster) + theme_grey(base_size = 8) + geom_point(aes(colour=df$cluster),size=3,alpha=0.3) + labs(title="K-means sobre componentes principales") + xlim(c(-5,5))  + ylim(c(-7,7)) + theme(panel.background = element_blank())


tables.cluster <- as.data.frame.matrix(table(df$cluster,df$position))
tables.cluster <- cbind(cluster=c(1,2,3,4),tables.cluster)

df.num$position <- df$position

plots.3d %<a-% scatterplot3d(df$pca1,df$pca2,df$pca3,angle=100,type="h",color=cc[df$position],pch=16,box=FALSE,xlab="standardized PC1",ylab="standardized PC2",zlab="standardized PC3")

# Análisis discriminante cuadrático

smp_size <- floor(0.75 * nrow(df.num))
train_ind <- sample(seq_len(nrow(df.num)), size = smp_size)
train <- df.num[train_ind, ]
test <- df.num[-train_ind, ]

model <- qda(position ~ ppv + bounce + time + traffic_mail + traffic_direct + traffic_display + traffic_social + traffic_social + traffic_referrals,data=train)

prediction <- predict(model,test)

test$score <- prediction$posterior[,2]
test$predicted <- prediction$class


tables.matrix <- confusionMatrix(test$predicted,test$position,positive="1",mode="sens_spec")
tables.matrix

?confusionMatrix

plots.qda <- ggplot(test,aes(x=score,colour=position,fill=position))+ theme_grey(base_size = 8) + geom_density(alpha=0.4) + theme(panel.background = element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +labs(title="Puntuaciones discriminantes")
