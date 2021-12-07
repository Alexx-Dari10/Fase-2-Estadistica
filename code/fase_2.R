library(lmtest)
library(rpart)

############################
##### Regresion Lineal #####
############################

### Datos
# variable dependiente: ResiduaryResistance
# variabes independientes: PrismCoeff, FroudeNumber 

yatch_data <- read.table('yatch.csv', header=TRUE, sep=',')
attach(yatch_data)

# analizando modelo regresion lineal multiple
multi.fit = lm( ResiduaryResistance ~ PrismCoeff + FroudeNumber, data = yatch_data )
summary(multi.fit)
## coeficiente prismatico ni el intercepto son significativos


# analizando matriz de correlacion
cor<-cor(yatch_data)
cor
## coeficiente prismatico no tiene correlacion con resistencia residual


# eliminando PrismCoeff del modelo
backward = lm( ResiduaryResistance ~ FroudeNumber, data = yatch_data )

# analizando nuevo modelo
summary(backward)


### Verificando supuestos

## Analizando los residuos
residuals <- backward$residuals

## 1. La media de los errores es cero y la suma de los errores es cero
mean(residuals)
sum(residuals)
# ok

## 2. Errores normalmente distribuidos
hist(residuals, main = "Histogram")

qqnorm(residuals)
qqline(residuals)

# test de Shapiro-Wilk
shapiro.test(residuals)
# no se cumple

## 3. Independencia de los residuos
dwtest(backward)
# no se cumple

## 4. Homocedasticidad
plot(backward$fitted.values, 
     rstandard(backward), 
     xlab = "Predictions", 
     ylab = "Residuals", 
     main = "Standardized Residuals")
abline(h=0,lty=2)

plot(ResiduaryResistance, 
     residuals, 
     xlab = "ResiduaryResistance", ylab = "Residuals", main = "ResiduaryResistance vs Residuals")
abline(h=0,lty=2)

bptest(backward)
# no se cumple




#################
##### ANOVA #####
#################

## datos
# variable dependiente: ResiduaryResistance
# factor: PrismCoeff
# bloques: FroudeNumber


# organizando datos
.data <- data.frame(PrismCoeff, FroudeNumber, ResiduaryResistance)

# comparando medias con graficos de caja y bigotes
boxplot(ResiduaryResistance ~ PrismCoeff, data = .data,
        main = "Boxplots. ResiduaryResistance vs PrismCoeff",
        ylim = c(0,30),
        xlab = "PrismCoeff",
        ylab = "ResiduaryResistance")

## medias bastante cercanas, es posible que el coeficiente prismatico no tenga efecto
## sobre la resistencia residual


boxplot(ResiduaryResistance ~ FroudeNumber, data = .data,
        main = "Boxplots. ResiduaryResistance vs FroudeNumber",
        xlab = "FroudeNumber",
        ylab = "ResiduaryResistance")

## es posible que el numnero de Froude tenga efecto sobre la Resistencia Residual


# analisis de ANOVA
anova <- aov(ResiduaryResistance ~ PrismCoeff + FroudeNumber, data = .data)
summary(anova)

## PrismCoeff p-valor > 0.05 => el coeficiente prismatico no influye en la resistencia 
# residual
## FroudeNumber p-valor << 0.05 => el numero de Froude influye en la resistencia
# residual


### Verificando supuestos

## Analizando los residuos
residuals <- anova$residuals

## 1. Errores normalmente distribuidos
hist(residuals, main = "Histogram")

qqnorm(residuals)
qqline(residuals)

# test de Shapiro-Wilk
shapiro.test(residuals)
# no se cumple

## 2. Independencia de los residuos
dwtest(anova)
# no se cumple

##3. Homocedasticidad
plot(anova$fitted.values, 
     rstandard(anova), 
     xlab = "Predictions", 
     ylab = "Residuals", 
     main = "Standardized Residuals")
abline(h=0,lty=2)

# prueba de bartlett
bartlett.test(residuals, FroudeNumber)
# no se cumple





##################################
##### Reduccion de dimension #####
##################################

### Analisis de componentes principales ###

# datos
data <- data.frame(LongPosition, PrismCoeff, LengthDispRatio, BeamDraughtRatio,
            LengthBeamRatio, FroudeNumber, ResiduaryResistance)

# analizando correlacion entre las variables
pairs(LongPosition ~ PrismCoeff + LengthDispRatio + BeamDraughtRatio +
        LengthBeamRatio + FroudeNumber + ResiduaryResistance, data = yatch_data, 
        main="Correlacion de las variables de la muestra")
## son demasiadas para analizar de esta forma

# matriz de correlacion en forma grafica
cor <- cor(data)
symnum(cor)
## no es una matriz altamente correlacionada

# analisis de componentes principales
acp <- prcomp(data, scale = TRUE)
summary(acp)
## nos podemos quedar con las 3 primeras componentes con un 71% aproxx

# graficando componentes principales 
plot(acp, main= "Plot de componentes principales")

# matriz de valores propios
acp$rotation

# biplot de las dos primeras componentes
biplot(acp)


### Cluster jerarquico ###

# estandarizando
data.std <- scale(data, scale = TRUE)

# estandarizados vs sin estandarizar
plot(data$LongPosition, data$PrismCoeff, main = "Datos sin estandarizar",
          xlab = "LongPosition", ylab = "PrismCoeff")
plot(data.std, main = "Datos estandarizados")

# matriz de distancia con la distancia euclidiana
d <- dist(data.std, method = "euclidean")

# ajuste completo
fit <- hclust(d, method = "complete")
d2 <- as.dendrogram(fit)

# graficando (2 clusers)
plot(fit)
rect.hclust(fit, k = 2, border = "red")

# graficando (4 clusters)
plot(fit)
rect.hclust(fit, k = 4, border = "red")


# por las medias
fit <- hclust(d, method = "average")
d2 <- as.dendrogram(fit)

# graficando (3 clusters)
plot(fit)
rect.hclust(fit, k = 3, border = "red")


### K-means ###

fit.k1 <- kmeans(data.std, 4)
fit.k1
## % de similitud de cada cluster bajo

fit.k2 <- kmeans(data.std, 8)
fit.k2
## 60 % no es lo ideal pero esta bien

# analisis grafico
pairs(data.std, col = fit.k2$cluster)

# analizando resistencia residual y coeficiente prismatico
plot(data.std[ ,'ResiduaryResistance'], data.std[ ,'PrismCoeff'], 
     col = fit.k2$cluster, lwd=2, ylab = "Residuary Resistance", xlab = "Prism Coeff")
points(fit.k2$centers, col = 1:12, pch = 6, lwd = 2)


### Arbol de clasificacion ###

# cardinalidad de la poblacion (308 pero puede variar)
l <- length(data[,1])
l

# sub = se escoge al azar las 2/3 partes de la poblacion
# como conjunto entrenante para crear el arbol
sub <- sample(1:l, 2*l/3)
sub

# conjunto entrenante
data[sub,]

# conjunto de prueba
data[-sub,]

# construyendo el arbol con el conjunto entrenante
# para predecir la resistencia residual a partir del resto de las variables
tree <- rpart(ResiduaryResistance ~ ., data = data[sub,], 
                  cp = 0.01, maxdepth = 3)

# sumario del arbol
summary(tree)
tree

# graficando arbol
plot(tree)
text(tree, use.n = TRUE, all = TRUE, pretty = 0, xpd = TRUE)

printcp(tree)
plotcp(tree)

# haciendo la prediccion con el conjunto de prueba
# teniendo en cuenta el arbol obtenido
pred <- predict(tree, newdata = data[-sub,], type = "vector")
pred

# matriz de confusion
tb <- table(pred, data[-sub,]$ResiduaryResistance)

# calculo del error
error <- 1 - (sum(diag(tb))/sum(tb))
error
