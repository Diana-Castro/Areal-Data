# TFM: ANÁLISIS DE DATOS AGREGADOS POR ÁREA- contaminación atmosférica y desigualdad socioeconómica
# Análisis de datos espaciales
# Unidad de análisis: 313 municipios de Galicia
# Contaminantes: NO2 y PM2.5 (datos CAMS/Copernicus) para el año 2023
# Datos socioeconómicos: Renta per cápita bruta en Galicia de: Instituto Galego de Estatística
# Autor: Diana Cristina Castro Armenta

# Librerías necesarias

library(terra)
library(exactextractr)
library(mapSpain)
library(sf)
library(dplyr)
library(tmap)
library(stringr)
library(ggplot2)
library(spdep)
library(leaflet)
library(patchwork) 
library(viridis)   
library(scales)
library(INLA)

# Leer datos
municipios<-st_read("data/galmun_contaminacion_renta.gpkg")
municipios_utm<- st_transform(municipios, crs= 32629) # Proyección UTM para análisis en metros
st_crs(municipios_utm)$units # Verificar CRS

# Estadísticas descriptivas
vars<- c("no2_median", "pm25_median", "renta.pc")
summary(st_drop_geometry(municipios[, vars]))

# 1. Construcción del vecindario y la matriz W
# 1.1 Designación vecindario
nb_queen<- poly2nb(municipios, queen = TRUE)# Vecindario

# Identificar municipios sin vecinos
cat("Municipios sin vecinos:", sum(card(nb_queen) == 0), "\n")
islas<- which(card(nb_queen) == 0)
cat("Municipios isla:", islas, "\n")
municipios$nombre_limpio[islas] #Correspondiente a la Isla de Arousa
coords<-st_centroid(st_geometry(municipios_utm)) |> st_coordinates()


# modificar para isla, manualmente, cambiar para hacer match automático 
# 1.1.1 Asignar el vecino más cercano por distancia a cada isla
nb_queen <- addlinks1(nb_queen,
                      from = islas,
                      to   = knn2nb(knearneigh(coords, k = 1))[[islas]])
# Agrega por orden de k: Cambados, Vilanova de Arousa, O Grove como vecinos
# municipios$name[nb_queen[[267]]] # A Estrada

# 1.1.2 Alternativa con criterio KNN para municipio aislado
nb_knn5  <-knn2nb(knearneigh(coords, k = 5))

# 1.2 Matriz de pesos w
Wq <- nb2listw(nb_queen, style ="W", zero.policy =TRUE) #?nb2listw

# 1.2.1 Matriz de pesos binaria no estandarizada (para modelos ICAR/BYM)
Wq_bin <- nb2listw(nb_queen, style= "B")

# 1.2.2 Matriz de pesos por inversa de la distancia (IDW) según KNN
distancias<- nbdists(nb_knn5, coords) #distancias entre vecinos
pesos_idw <-lapply(distancias, function(d) 1/d) #pesos por inversa de la distancia
# Lista de pesos estandarizada por filas
Widw <- nb2listw(nb_knn5,
                 glist= pesos_idw,
                 style ="W")  # estandarización por filas
# 1.2.3 Verificar pesos
print(Wq)
summary(Wq)

print(Widw)

# 1.3 Visualización del grafo del vecindario
windows()
par(mfrow= c(1, 2))
plot(st_geometry(municipios_utm), border = "grey70",main ="Vecindad Queen")
plot(nb_queen, coords,add =TRUE, col ="steelblue",lwd =0.4)

# Vecinos de ese municipio (índices reales dentro del nb completo)
vecinos_267 <- nb_queen[[267]]
# Segmentos manuales municipio -> cada vecino
for (v in vecinos_267) {
  segments(coords[267, 1], coords[267, 2],
           coords[v, 1], coords[v, 2],
           col = "red", lwd = 2)
}
# resaltar A Estrada
plot(st_geometry(municipios_utm), border ="grey70", main= "Vecindad KNN (k=5)")
plot(nb_knn5, coords, add= TRUE, col= "tomato",lwd = 0.4)
par(mfrow=c(1, 1))
# Se utilizará matriz queen ajustada para análisis de autocorrelación espacial y 
#modelos espaciales, ya que refleja mejor la contigüidad entre municipios y el 
#tipo de fenómeno que es la dispersión de contaminates atmosféricos.


# 2. Análisis exploratorio de datos espaciales (ESDA)
# Paleta de colores
paleta <- c('#5e3c99', '#b2abd2', '#fdb863', '#e66101')

# Paleta para gradientes continuos (de bajo a alto)
paleta_continua_low  <- '#5e3c99'   # morado oscuro = valores bajos
paleta_continua_high <- '#e66101'   # naranja oscuro = valores altos

# Paleta para LISA y BiLISA (4 cuadrantes + No significativo)
# Alto-Alto = naranja oscuro (zona crítica)
# Bajo-Bajo = morado oscuro (zona favorable)
# Alto-Bajo = naranja claro (outlier)
# Bajo-Alto = morado claro (outlier)
# No significativo = gris claro 
paleta_lisas <- c(
  "Alto-Alto"        = '#e66101',
  "Bajo-Bajo"        = '#5e3c99',
  "Alto-Bajo"        = '#fdb863',
  "Bajo-Alto"        = '#b2abd2',
  "No significativo" = 'grey97'
)
# 2.1 Mapas coropléticos
# mostrando solo variable continua sin discretizar
library(ggplot2)

p_no2<- ggplot(municipios) + geom_sf(aes(fill= no2_median), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient(low = paleta_continua_low, high = paleta_continua_high,
                      name = "NO2 (µg/m3)") + 
  theme_minimal() + labs(title = "NO2 mediana (µg/m3)" ) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
p_pm25<- ggplot(municipios) + geom_sf(aes(fill =pm25_median), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient(low = paleta_continua_low, high = paleta_continua_high,
                      name = "PM2.5 (µg/m3)") +
  theme_minimal() + labs(title = "PM2.5 mediana (µg/m3)") +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
p_renta<- ggplot(municipios) + geom_sf(aes(fill=renta.pc), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient(low = paleta_continua_low, high = paleta_continua_high,
                      name = "Renta (miles €/cápita)") +
  theme_minimal() + labs(title ="Renta bruta (miles de euros per cápita)") +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
windows()
p_pm25 | p_no2 |p_renta +
  plot_annotation(
    title = "Distribución espacial de contaminantes y renta en Galicia", 
  )

# Boxplots

bp_pm25<- ggplot(municipios, aes(x = "", y = pm25_median)) +
  geom_boxplot(fill = paleta_continua_high, color = "grey30") +
  labs(title = "Boxplot de PM2.5 mediana", y = "PM2.5 (µg/m3)", x = "") +
  theme_minimal()
bp_no2<- ggplot(municipios, aes(x = "", y = no2_median)) +
  geom_boxplot(fill = paleta_continua_high, color = "grey30") +
  labs(title = "Boxplot de NO2 mediana", y = "NO2 (µg/m3)", x = "") +
  theme_minimal()
bp_renta<- ggplot(municipios, aes(x = "", y = renta.pc)) +
  geom_boxplot(fill = paleta_continua_high, color = "grey30") +
  labs(title = "Boxplot de Renta bruta per cápita", y = "Renta (miles €/cápita)", x = "") +
  theme_minimal()
windows()
bp_pm25 | bp_no2 | bp_renta +
  plot_annotation(
    title = "Boxplots de contaminantes y renta en Galicia",
  )

# 2.2 Autocorrelación espacial global: Índice I de Moran
EI<--1/(313-1) #esperanza de I Moran 
semilla<- 27750018 # para reproducibilidad
set.seed(semilla)
# PM2.5
moran_pm25<-moran.test(municipios$pm25_median,
                       listw=Wq,
                       alternative ="two.sided") # para detectar cualquier tipo de autocorrelación espacial
moran_pm25_perm <- moran.mc(municipios$pm25_median,
                            listw =Wq,
                            nsim=999, alternative ="two.sided")
#moran.test(municipios$pm25_median, listw = Widw, alternative = "two.sided")
# NO2
moran_no2<-moran.test(municipios$no2_median,
                      listw= Wq,
                      alternative="two.sided")
moran_no2_perm <- moran.mc(municipios$no2_median,
                            listw=Wq,
                            nsim= 999, alternative="two.sided")
# Renta
moran_renta<-moran.test(municipios$renta_log,
                         listw = Wq,
                         alternative ="two.sided")
moran_renta_perm <- moran.mc(municipios$renta_log,
                             listw= Wq,
                             nsim= 999, alternative="two.sided")
# Resultados
cat("\n--- Índice de Moran global (Queen, estand. filas) ---\n")
cat(sprintf("PM2.5: I = %.4f, p-valor= %.4f\n",
            moran_pm25$estimate[1], moran_pm25_perm$p.value))
cat(sprintf("NO2: I = %.4f, p-valor= %.4f\n",
            moran_no2$estimate[1],  moran_no2_perm$p.value))
cat(sprintf("Log-renta: I = %.4f, p-valor= %.4f\n",
            moran_renta$estimate[1], moran_renta_perm$p.value))
# Para todas las variables se rechaza la hipótesis nula de no autocorrelación espacial, 
# indicando que hay una estructura espacial global significativa.

# 2.3 Diagrama de dispersión de Moran Global
windows()
par(mfrow= c(1, 2))

# PM2.5
moran.plot(scale(municipios$pm25_median, scale=FALSE)[,1], #Evitar escalar, solo centrar para diagrama de Moran
           listw = Wq,
           labels= municipios$name,
           xlab= "PM2.5 (centrado)",
           ylab="Retardo espacial de PM2.5",
           main="Diagrama de dispersión de Moran — PM2.5")
# NO2

moran.plot(scale(municipios$no2_median, scale=FALSE)[,1],
           listw = Wq,
           labels= municipios$name,
           xlab= "NO2 (centrado)",
           ylab= "Retardo espacial de NO2",
           main="Diagrama de dispersión de Moran — NO2")
# Renta
moran.plot(scale(municipios$renta_log, scale=FALSE)[,1],
           listw = Wq,
           labels= municipios$name,
           xlab="Log-renta (centrado)",
           ylab= "Retardo espacial de log-renta",
           main= "Diagrama de dispersión de Moran — Log-renta")

# 2.4 Autocorrelación espacial local: Índice I Local de Moran 
set.seed(semilla)
# PM2.5
locm_pm25<-localmoran(x=municipios$pm25_median, listw = Wq) # alternative= "two.sided" por default
head(locm_pm25) # dist.normal
locm_pm25_perm<-localmoran_perm(x=municipios$pm25_median, listw = Wq, nsim=999)
head(locm_pm25_perm) # dist. de permuaciones
# NO2
locm_no2<-localmoran(x=municipios$no2_median, listw = Wq)
locm_no2_perm<-localmoran_perm(x=municipios$no2_median, listw = Wq, nsim=999)
# Renta
locm_renta<-localmoran(x=municipios$renta_log, listw =Wq)
locm_renta_perm<-localmoran_perm(x=municipios$renta_log, listw=Wq, nsim=999)
# --2.4.1 Agregar el índice de Moran local al objeto municipios
municipios$I_localpm25 <- locm_pm25[,"Ii"]
#summary(locm_pm25[,"Ii"])
# Variable categorizada
municipios$I_local_catpm25 <- cut(
  municipios$I_localpm25,
  breaks = c(-Inf, 0.08, 0.5, 1.3, Inf),
  labels = c(
    "Negativo alto",
    "Negativo bajo",
    "Positivo bajo",
    "Positivo alto"
  ),
  right = FALSE
)
municipios$I_lmz_pm25 <- locm_pm25[,"Z.Ii"]
municipios$lmp_pm25 <- locm_pm25[,"Pr(z != E(Ii))"] # p-valor 

# 2.4.2 Representación del índice de Moran local para PM2.5
p1 <- ggplot(municipios) +
  geom_sf(aes(fill = I_local_catpm25), colour = "grey40", linewidth = 0.15) +
  scale_fill_manual(values = paleta,
                    name = "Índice local\nde Moran",
                    na.value = "grey92") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
p2 <- ggplot(municipios) +
  geom_sf(aes(fill = pm25_median), colour = "grey40", linewidth = 0.15) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey97",
                      high = paleta_continua_high, midpoint = median(municipios$pm25_median, na.rm = TRUE),
                      name = "Concentración PM 2.5") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
p3 <- ggplot(municipios) +
  # incluir breaks para categorizar en qnorm(0.975) = 1.96 y qnorm(0.025) = -1.96
  geom_sf(aes(fill = I_lmz_pm25), colour = "grey40", linewidth = 0.15) +
  scale_x_continuous(breaks = c(-Inf, -1.96, 1.96, Inf)) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey97",
                       high = paleta_continua_high, midpoint = 0,
                       name = "Z-score\nÍndice local Moran") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

p4<- ggplot(municipios) +
  geom_sf(aes(fill = lmp_pm25), colour = "grey40", linewidth = 0.15) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey97",
                      high = paleta_continua_high, midpoint = 0.5,
                      name = "p-valor\n") +
  scale_x_continuous(breaks = c(-Inf, 0.05, Inf)) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
  )
windows()
(p1|p2) / (p3|p4) +
  plot_annotation(
    title = "Índice local de Moran para PM 2.5 en Galicia",
    subtitle = "313 municipios"
  )

# -- 2.4.3 Representación de clústers 
# Función auxiliar: clasificación en cuadrantes
clasificar_cuadrantes <- function(z, lw, lisa, alpha = 0.05) {
  z_c   <- z - mean(z, na.rm = TRUE) #variable centrada
  lag_z <- lag.listw(x=lw, var=z_c)
  quad <- rep(NA_character_, length(z))
  quad  <- case_when(
    z_c > 0 & lag_z > 0 ~ "Alto-Alto",
    z_c < 0 & lag_z < 0 ~ "Bajo-Bajo",
    z_c > 0 & lag_z < 0 ~ "Alto-Bajo",
    z_c < 0 & lag_z > 0 ~ "Bajo-Alto"
  )
  # Solo áreas significativas
  quad_sig <- quad
  quad_sig[lisa[, "Pr(z != E(Ii)) Sim"] > alpha] <- NA_character_ # Sim porque es para permutaciones
  list(quad = quad, quad_sig = quad_sig)
}

# Centrar variables
pm25_c <- municipios$pm25_median - mean(municipios$pm25_median)
no2_c  <- municipios$no2_median  - mean(municipios$no2_median)
rental_c <- municipios$renta_log - mean(municipios$renta_log)
# Clasificar según cuadrantes
res_pm25 <- clasificar_cuadrantes(pm25_c, Wq, locm_pm25_perm) # para permutaciones
res_no2  <- clasificar_cuadrantes(no2_c,  Wq, locm_no2_perm)
res_renta <- clasificar_cuadrantes(rental_c, Wq, locm_renta_perm)

niveles <- c("Alto-Alto", "Bajo-Bajo", "Alto-Bajo", "Bajo-Alto")


municipios <- municipios |>
  mutate(
    LISA_PM25     = factor(res_pm25$quad,     levels = niveles),
    LISA_PM25_SIG = factor(res_pm25$quad_sig, levels = niveles),
    LISA_NO2      = factor(res_no2$quad,      levels = niveles),
    LISA_NO2_SIG  = factor(res_no2$quad_sig,  levels = niveles),
    LISA_RENTA    = factor(res_renta$quad,    levels = niveles),
    LISA_RENTA_SIG= factor(res_renta$quad_sig,levels = niveles)
  )
# Mapas con clasificación LISA 
p_lisa_pm25 <- ggplot(municipios) +
  geom_sf(aes(fill = LISA_PM25), colour = "grey24", linewidth = 0.05) +
  scale_fill_manual(values = paleta_lisas, na.value = "grey90",
                    name = "Clasificación LISA", na.translate = TRUE, drop = FALSE) +
  labs(title = "LISA PM2.5") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

p_lisa_no2 <- ggplot(municipios) +
  geom_sf(aes(fill = LISA_NO2), colour = "grey24", linewidth = 0.05) +
  scale_fill_manual(values = paleta_lisas, na.value = "grey90",
                    name = "Clasificación LISA", na.translate = TRUE, drop = FALSE) +
  labs(title = "LISA NO2 ") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

p_lisa_renta <- ggplot(municipios) +
  geom_sf(aes(fill = LISA_RENTA), colour = "grey24", linewidth = 0.05) +
  scale_fill_manual(values = paleta_lisas, na.value = "grey90",
                    name = "Clasificación LISA", na.translate = TRUE, drop = FALSE) +
  labs(title = "LISA Log(Renta)") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

windows()
p_lisa_pm25 | p_lisa_no2 | p_lisa_renta +
  plot_annotation(
    title = "Clústeres LISA significativos para PM2.5, NO2 y Renta en Galicia",
    subtitle = "313 municipios"
  )

# Mapas LISA significativos
p_lisasig_pm25 <- ggplot(municipios) +
  geom_sf(aes(fill = LISA_PM25_SIG), colour = "grey24", linewidth = 0.05) +
  scale_fill_manual(values = paleta_lisas, na.value = "grey90",
                    name = "Clasificación (sig.)", na.translate = TRUE, drop = FALSE,
                    labels = c(niveles, "No significativo")) +
  labs(title = "LISA PM2.5 (p < 0.05, permutaciones)") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

p_lisasig_no2 <- ggplot(municipios) +
  geom_sf(aes(fill = LISA_NO2_SIG), colour = "grey24", linewidth = 0.05) +
  scale_fill_manual(values = paleta_lisas, na.value = "grey90",
                    name = "Clasificación (sig.)", na.translate = TRUE, drop = FALSE,
                    labels = c(niveles, "No significativo")) +
  labs(title = "LISA NO2 (p < 0.05, permutaciones)") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
p_lisasig_renta <- ggplot(municipios) +
  geom_sf(aes(fill = LISA_RENTA_SIG), colour = "grey24", linewidth = 0.05) +
  scale_fill_manual(values = paleta_lisas, na.value = "grey90",
                    name = "Clasificación (sig.)", na.translate = TRUE, drop = FALSE,
                    labels = c(niveles, "No significativo")) +
  labs(title = "LISA log(Renta) (p < 0.05, permutaciones)") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
windows()
p_lisasig_pm25 | p_lisasig_no2 | p_lisasig_renta +
  plot_annotation(
    title = "Clústeres LISA significativos para PM2.5, NO2 y Renta en Galicia",
    subtitle = "313 municipios"
  )
#ggsave("mapassign_lisa.png", width = 14, height = 6, dpi = 300)

# Tabla de frecuencias
(tabla_lisa_pm25 <- table(municipios$LISA_PM25_SIG))
(tabla_lisa_no2 <- table(municipios$LISA_NO2_SIG))
(tabla_lisa_renta <- table(municipios$LISA_RENTA_SIG))
which(municipios$LISA_RENTA_SIG=="Alto-Bajo")

#########################################################
# 3.5 BiLISA
library(bispdep)
# Estandarizar variables para modelos lineales
municipios$pm25_std <- as.numeric(scale(municipios$pm25_median)[,1])
municipios$no2_std  <- as.numeric(scale(municipios$no2_median)[,1])
municipios$logrenta_std <- as.numeric(scale(municipios$renta_log)[,1])

# Retardos espaciales
lag_logrenta <- lag.listw(Wq, municipios$logrenta_std)
lag_pm25 <- lag.listw(Wq, municipios$pm25_std)
lag_no2 <- lag.listw(Wq, municipios$no2_std)

# Índice Moran Global
(moranbi_pm25 <- moran.bi(municipios$pm25_std, lag_logrenta, listw = Wq))
(moranbi_no2 <- moran.bi(municipios$no2_std, lag_logrenta, listw = Wq))


# Diagrama de dispersión de Moran bivariante global
#moranbi.plot(varX = municipios$pm25_std, varY = municipios$logrenta_std, listw = Wq,
#             xlab = "PM2.5 (estandarizado)", ylab = "Retardo espacial de log-renta",
#            main = "Diagrama de dispersión BiLISA: PM2.5 x Renta bruta per capita")

#moranbi.plot(varX = municipios$no2_std, varY = municipios$logrenta_std, listw = Wq,
#             xlab = "NO2 (estandarizado)", ylab = "Retardo espacial de log-renta",
#             main = "Diagrama de dispersión BiLISA: NO2 x Renta bruta per capita")


# BiLISAs 
bilisa_pm25 <- localmoran.bi(varX = municipios$pm25_std, varY = municipios$logrenta_std, listw = Wq)
# alternative = "two.sided" por default
bilisa_no2 <- localmoran.bi(varX = municipios$no2_std, varY = municipios$logrenta_std, listw = Wq)

# Mapas de Clústeres BiLISA- DON'T WORK
moranbi.cluster(varY = municipios$pm25_std, varX = municipios$logrenta_std, listw = Wq, 
                # significant = TRUE y alternative = "two.sided" por default
                pleg = "topleft", polygons = municipios$geom)
title("BiLISA: PM2.5 x Renta bruta per capita - Galicia")
moranbi.cluster(varY = municipios$no2_std, varX = municipios$logrenta_std, listw = Wq, 
                # significant = TRUE y alternative = "two.sided" por default
                pleg = "topleft", polygons = municipios$geom)
title("BiLISA: NO2 x Renta bruta per capita - Galicia")


# Municipios con significativamente alta contaminación y renta baja (High-Low)
# seguidos de municipios con baja contaminación y renta alta (Low-High)

municipios$name[which(municipios$pm25_std > 0 & lag_logrenta < 0 & bilisa_pm25[,"Pr(z != 0)"] < 0.05)]
municipios$name[which(municipios$pm25_std < 0 & lag_logrenta > 0 & bilisa_pm25[,"Pr(z != 0)"] < 0.05)]

municipios$name[which(municipios$no2_std > 0 & lag_logrenta < 0 & bilisa_no2[,"Pr(z != 0)"] < 0.05)]
municipios$name[which(municipios$no2_std < 0 & lag_logrenta > 0 & bilisa_no2[,"Pr(z != 0)"] < 0.05)]
# no coinciden las clasificaciones de los municipioscon moranbi.cluster
# Manualmente

municipios <- municipios |>
  mutate(
    bilisa_pm25_quad = case_when(
      pm25_std > 0 & lag_logrenta < 0 & bilisa_pm25[, "Pr(z != 0)"] < 0.05 ~ "Alto-Bajo",
      pm25_std < 0 & lag_logrenta > 0 & bilisa_pm25[, "Pr(z != 0)"] < 0.05 ~ "Bajo-Alto",
      pm25_std > 0 & lag_logrenta > 0 & bilisa_pm25[, "Pr(z != 0)"] < 0.05 ~ "Alto-Alto",
      pm25_std < 0 & lag_logrenta < 0 & bilisa_pm25[, "Pr(z != 0)"] < 0.05 ~ "Bajo-Bajo",
      TRUE ~ "No significativo"
    ),
    bilisa_no2_quad = case_when(
      no2_std > 0 & lag_logrenta < 0 & bilisa_no2[, "Pr(z != 0)"] < 0.05 ~ "Alto-Bajo",
      no2_std < 0 & lag_logrenta > 0 & bilisa_no2[, "Pr(z != 0)"] < 0.05 ~ "Bajo-Alto",
      no2_std > 0 & lag_logrenta > 0 & bilisa_no2[, "Pr(z != 0)"] < 0.05 ~ "Alto-Alto",
      no2_std < 0 & lag_logrenta < 0 & bilisa_no2[, "Pr(z != 0)"] < 0.05 ~ "Bajo-Bajo",
      TRUE ~ "No significativo"
    )
  )

niveles_bilisa <- c("Alto-Alto","Bajo-Bajo", "Alto-Bajo", "Bajo-Alto", "No significativo")

# Tabla de frecuencias
(tabla_bilisa_pm25 <- table(municipios$bilisa_pm25_quad))
(tabla_bilisa_no2 <- table(municipios$bilisa_no2_quad))
which(municipios$bilisa_pm25_quad=="Alto-Bajo")

# Tabla municipios BiLISA Alto-Bajo 
tabla_bilisa_ab <- municipios |>
  st_drop_geometry() |>
  filter(bilisa_pm25_quad == "Alto-Bajo" | bilisa_no2_quad == "Alto-Bajo") |>
  transmute(
    Municipio = name,
    Provincia = ine.prov.name,  
    PM2.5_AltoBajo = bilisa_pm25_quad == "Alto-Bajo",
    NO2_AltoBajo   = bilisa_no2_quad == "Alto-Bajo"
  )
print(tabla_bilisa_ab)


# Exportar directamente a LaTeX
#library(xtable)
#print(xtable(tabla_bilisa_ab, caption = "Municipios en el cuadrante BiLISA Alto-Bajo"),
#      include.rownames = FALSE)
# Mapa PM2.5
p_bilisa_pm25 <- ggplot(municipios) +
  geom_sf(aes(fill = factor(bilisa_pm25_quad, levels = niveles_bilisa)),
          colour = "grey24", linewidth = 0.05) +
  scale_fill_manual(values = paleta_lisas, name = "BiLISA",
                    na.value = "grey99") +
  labs(title    = "BiLISA: PM2.5 x Log(Renta) — Galicia",
       subtitle = "p < 0.05 (permutaciones)") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

# Mapa NO2
p_bilisa_no2 <- ggplot(municipios) +
  geom_sf(aes(fill = factor(bilisa_no2_quad, levels = niveles_bilisa)),
          colour = "grey24", linewidth = 0.05) +
  scale_fill_manual(values = paleta_lisas, name = "BiLISA",
                    na.value = "grey99") +
  labs(title    = "BiLISA: NO2 x Log(Renta) — Galicia",
       subtitle = "p < 0.05 (permutaciones)") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

windows()
p_bilisa_pm25 | p_bilisa_no2


#########################################################



# --- 3. Modelización
library(GGally)
windows()
ggpairs(data=municipios, columns = c("pm25_median", "no2_median", "renta_log"))
# Extraer la matriz de correlación numérica 
cor(st_drop_geometry(municipios)[, c("pm25_median", "no2_median", "renta_log")],
    use = "complete.obs")

# 3.1 Modelos no espaciales (baseline)
# Modelos lineales ordinarios: contaminantes en función de log-renta 
lm_pm25 <- lm(pm25_median ~ renta_log, data = municipios) 
summary(lm_pm25)

windows()
plot(lm_pm25)
lm_no2  <- lm(no2_median ~ renta_log, data = municipios)
summary(lm_no2)
# en ambos modelos p-valor para coef. de renta es significativo y positivo,
# indicando que a mayor renta per cápita, mayor concentración de contaminantes.
# En el caso del NO2 tanto intercepto como coeficiente de renta_log es significativo

#3.1.1 Diagnóstico de autocorrelación en residuos (H_0:No autocorrelación)

municipios$resid_pm25 <- residuals(lm_pm25)
municipios$resid_no2  <- residuals(lm_no2)
set.seed(semilla)
moran_res_pm25 <- moran.mc(municipios$resid_pm25, listw = Wq,
                           alternative = "two.sided", nsim = 999)
moran_res_no2<- moran.mc(municipios$resid_no2, listw = Wq,
                           alternative = "two.sided", nsim = 999)
cat("\n--- Moran de residuos ---\n")
cat(sprintf("PM2.5: I = %.4f, p-valor = %.4f\n",
            moran_res_pm25$statistic, moran_res_pm25$p.value))
cat(sprintf("NO2  : I = %.4f, p-valor = %.4f\n",
            moran_res_no2$statistic,  moran_res_no2$p.value))
#Se rechaza H_0, indicando presencia de autocorrelación en residuos

# representación en mapa de residuos
p_res_pm25 <- ggplot(municipios) +
  geom_sf(aes(fill = resid_pm25), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey97",
                       high = paleta_continua_high, midpoint = 0,
                       name = "Residuos PM2.5") +
  labs(title = "Residuos del modelo lineal para PM2.5") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
p_res_no2 <- ggplot(municipios) +
  geom_sf(aes(fill = resid_no2), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey92",
                       high = paleta_continua_high, midpoint = 0,
                       name = "Residuos PM2.5") +
  labs(title = "Residuos del modelo lineal para NO2") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
windows()
p_res_pm25 | p_res_no2


#ggsave("outputs/residuos_ols.png", width = 14, height = 6, dpi = 300)

# Valores de I cercanos a 1, indicando todavía fuerte autocorrelación espacial positiva
# o presencia de dependencia espacial residual



# 3.2.1 Modelo Gaussiano de referencia con INLA
lm_pm25_inla <- inla(formula= pm25_median ~ renta_log , 
                     data = municipios,
                     family ="gaussian",
                     control.compute = list(dic = TRUE, waic = TRUE),
                     control.predictor = list(compute = TRUE))
summary(lm_pm25_inla)
lm_no2_inla <- inla(formula= no2_median ~ renta_log, 
                    data = municipios,
                    family ="gaussian",
                    control.compute = list(dic = TRUE, waic = TRUE),
                    control.predictor = list(compute = TRUE))
summary(lm_no2_inla)
# En ambos modelos, el coeficiente de log-renta coincide con el generado por lm.

#3.2.2. Modelo INLA tipo BYM (Bayesiano)
#   Z_i = beta_0 + beta_1 * LOG_RENTA_STD_i + U_i + V_i
#   U ~ ICAR según W    [comp. esp. estructurada]
#   V_i ~ N(0, sigma2_v)    [comp. esp. no estructurada]
# Matriz para INLA no estandarizada tipo queen
Wq_bin
W_inla <- as(Wq_bin, "CsparseMatrix") # formato para INLA

# Indentificadores de área
municipios$id <- 1:nrow(municipios)
# Modelo espacial para PM2.5
sp_pm25_inla <- inla(formula = pm25_median ~ renta_log + f(id, model = "besagproper", graph = W_inla),
                     data = st_drop_geometry(municipios),
                     family = "gaussian",
                     control.compute = list(dic = TRUE, waic = TRUE , return.marginals.predictor = TRUE),
                     control.predictor = list(compute = TRUE))
summary(sp_pm25_inla)
# disminuye el coeficiente de log-renta, de 1.153 a 0.136, pero sigue siendo significativo: IC no incluye cero
# Modelo espacial para NO2
sp_no2_inla <- inla(formula = no2_median ~ renta_log + f(id, model = "besagproper", graph = W_inla),
                    data = st_drop_geometry(municipios),
                    family = "gaussian",
                    control.compute = list(dic = TRUE, waic = TRUE, return.marginals.predictor = TRUE),
                    control.predictor = list(compute = TRUE))
summary(sp_no2_inla)
# coeficiente de log_renta=0.384, antes 2.884, caso similar a pm2.5, pero ahora intercept ya no es significativo

# Proporción de varianza explicada por componente espacial (U+V) frente a la varianza total
# con la proporción de varianza espacial estructurada (phi)
# Nombres exactos de los hiperparámetros del modelo:
names(sp_pm25_inla$marginals.hyperpar)

cal_phi <- function(modelo_inla) {
  marg_prec_U <- modelo_inla$marginals.hyperpar$`Precision for id`
  marg_prec_V <- modelo_inla$marginals.hyperpar$`Precision for the Gaussian observations`
  
  var_U <- inla.emarginal(function(x) 1/x, marg_prec_U)
  var_V <- inla.emarginal(function(x) 1/x, marg_prec_V)
  phi   <- var_U / (var_U + var_V)
  
  # IC 95% de phi vía muestreo de las marginales
  set.seed(semilla)
  sim_prec_U <- inla.rmarginal(10000, marg_prec_U)
  sim_prec_V <- inla.rmarginal(10000, marg_prec_V)
  sim_phi <- (1/sim_prec_U) / (1/sim_prec_U + 1/sim_prec_V)
  
  list(var_U = var_U, var_V = var_V, phi = phi,
       phi_ic = quantile(sim_phi, c(0.025, 0.975)))
}

phi_pm25 <- cal_phi(sp_pm25_inla)
phi_no2  <- cal_phi(sp_no2_inla)


cat(sprintf("PM2.5: phi = %.4f (IC95%%: %.4f, %.4f)\n",
            phi_pm25$phi, phi_pm25$phi_ic[1], phi_pm25$phi_ic[2]))
cat(sprintf("NO2  : phi = %.4f (IC95%%: %.4f, %.4f)\n",
            phi_no2$phi, phi_no2$phi_ic[1], phi_no2$phi_ic[2]))


# --3.2.3 Extracción de componentes
# Componentes espaciales estructuradas (U)
municipios$u_pm25 <- sp_pm25_inla$summary.random$id$mean
range(municipios$u_pm25)
municipios$u_no2  <- sp_no2_inla$summary.random$id$mean
range(municipios$u_no2)

# Componentes espaciales no estructuradas (V)
# V= observados -predictor lineal (incluye U)
municipios$v_pm25 <- municipios$pm25_median - sp_pm25_inla$summary.linear.predictor$mean
range(municipios$v_pm25)
municipios$v_no2  <- municipios$no2_median - sp_no2_inla$summary.linear.predictor$mean
range(municipios$v_no2)  

# Representación coroplética
# Componentes espaciales estructuradas (U)
p_u_pm25 <- ggplot(municipios) +
  geom_sf(aes(fill = u_pm25), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey97",
                       high = paleta_continua_high, midpoint = 0,
                       name = "Componente estructurada (U)",
                       limits = range(municipios$u_pm25)) +
  labs(title = "Componente espacial estructurada para PM2.5") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
p_u_no2 <- ggplot(municipios) +
  geom_sf(aes(fill = u_no2), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey97",
                       high = paleta_continua_high, midpoint = 0,
                       name = "Componente estructurada (U)",
                       limits = range(municipios$u_no2)) +
  labs(title = "Componente espacial estructurada para NO2") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
# Componentes no estructuradas (V)
p_v_pm25 <- ggplot(municipios) +
  geom_sf(aes(fill = v_pm25), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey97",
                       high = paleta_continua_high, midpoint = 0,
                       name = "Componente no estructurada (V)",
                       limits = range(municipios$v_pm25)) +
  labs(title = "Componente espacial no estructurada para PM2.5") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
p_v_no2 <- ggplot(municipios) +
  geom_sf(aes(fill = v_no2), colour = "grey24", linewidth = 0.05) +
  scale_fill_gradient2(low = paleta_continua_low, mid = "grey97",
                       high = paleta_continua_high, midpoint = 0,
                       name = "Componente no estructurada (V)",
                       limits = range(municipios$v_no2)) +
  labs(title = "Componente espacial no estructurada para NO2") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
windows()
(p_u_pm25 | p_u_no2) / (p_v_pm25 | p_v_no2) +
  plot_annotation(
    title = "Componentes espaciales del modelo BYM para PM2.5 y NO2 en Galicia",
    subtitle = "313 municipios")

# TEST de Moran MC
set.seed(semilla)
moran.mc(x=municipios$v_pm25, listw=Wq, alternative = "two.sided", nsim = 999)
moran.mc(x=municipios$v_no2, listw=Wq, alternative = "two.sided", nsim = 999)

# Ambos contaminantes muestran disminución en correlación espacial y visualmente la componente no estructurada 
# pierde patrón espacial, p-valor no significativo :. no se rechaza H0: no autocorrelación espacial, indicando que el modelo BYM ha capturado la dependencia espacial en los datos.
  
# --- 3.3. Análisis de sensibilidad de vecindario KNN (k=5) vs queen
nb_knn5_sym <- make.sym.nb(nb_knn5) # es necesario simetrizar la matriz de pesos para INLA, por lo que se fuerza simetría
Wknn_bin <- nb2listw(nb_knn5_sym, style = "B")
W_knn_inla <- as(Wknn_bin, "CsparseMatrix")

sp_pm25_knn <- inla(formula = pm25_median ~ renta_log + f(id, model = "besagproper", graph = W_knn_inla),
                    data = st_drop_geometry(municipios),
                    family = "gaussian",
                    control.compute = list(dic = TRUE, waic = TRUE),
                    control.predictor = list(compute = TRUE))

sp_no2_knn <- inla(formula = no2_median ~ renta_log + f(id, model = "besagproper", graph = W_knn_inla),
                   data = st_drop_geometry(municipios),
                   family = "gaussian",
                   control.compute = list(dic = TRUE, waic = TRUE),
                   control.predictor = list(compute = TRUE))

# Tabla comparativa de beta_1 (log-renta) bajo ambos criterios de vecindario
tabla_sensibilidad <- data.frame(
  Modelo = c("PM2.5", "NO2"),
  beta_queen = c(sp_pm25_inla$summary.fixed["renta_log", "mean"],
                 sp_no2_inla$summary.fixed["renta_log", "mean"]),
  beta_queen_lo = c(sp_pm25_inla$summary.fixed["renta_log", "0.025quant"],
                    sp_no2_inla$summary.fixed["renta_log", "0.025quant"]),
  beta_queen_hi = c(sp_pm25_inla$summary.fixed["renta_log", "0.975quant"],
                    sp_no2_inla$summary.fixed["renta_log", "0.975quant"]),
  beta_knn = c(sp_pm25_knn$summary.fixed["renta_log", "mean"],
               sp_no2_knn$summary.fixed["renta_log", "mean"]),
  beta_knn_lo = c(sp_pm25_knn$summary.fixed["renta_log", "0.025quant"],
                  sp_no2_knn$summary.fixed["renta_log", "0.025quant"]),
  beta_knn_hi = c(sp_pm25_knn$summary.fixed["renta_log", "0.975quant"],
                  sp_no2_knn$summary.fixed["renta_log", "0.975quant"])
)
print(tabla_sensibilidad)
  
# 3.4 Evaluación de excedencia a directrices de la OMS
umbral_pm25 <- 5   # µg/m3, directriz OMS
umbral_no2  <- 10  # µg/m3, directriz OMS

prob_exceed_pm25 <- sapply(sp_pm25_inla$marginals.fitted.values, function(m) {
  1 - inla.pmarginal(umbral_pm25, m)
})
municipios$prob_exceed_no2 <- sapply(sp_no2_inla$marginals.fitted.values, function(m) {
  1 - inla.pmarginal(umbral_no2, m)
})

summary(municipios$prob_exceed_pm25)
summary(municipios$prob_exceed_no2)
# Municipios con probabilidad de excedencia relevante (>10%, ajustar umbral según resultado)
municipios$name[municipios$prob_exceed_pm25 > 0.10]
municipios$name[municipios$prob_exceed_no2  > 0.10]
