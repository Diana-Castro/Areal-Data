# ANÁLISIS ESPACIAL: CONTAMINACIÓN ATMOSFÉRICA Y DESIGUALDAD SOCIOECONÓMICA
# Unidad de análisis: 313 municipios de Galicia
# Contaminantes: NO2 y PM2.5 (datos CAMS/Copernicus) para el año 2023
# Datos socioeconómicos: Renta primaria bruta en Galicia de: Instituto Galego de Estatística
# Autor: Diana Cristina Castro Armenta

# Librerías necesarias

library(terra)
library(exactextractr)
library(mapSpain)
library(sf)
library(dplyr)

# 1. Límites Administrativos: municipios de Galicia
municipios <- st_transform(esp_get_munic(region = "ES11"), crs= 4326) # municipios_galicia <- esp_get_munic() |> filter(codauto == "12") |> st_transform(crs = 4326)
# Verificar
print(municipios)
windows()
plot(st_geometry(municipios), main = "313 Municipios de Galicia")


# 2. Cargar el raster NetCDF de datos CAMS
meses <- sprintf("%02d", 1:12) 
contaminantes <- c("no2", "pm2p5")

# Función que calcula la media mensual
media_mensual <- function(mes, cont) {
  ruta <- file.path("datos", 
                    paste0("cams.eaq.vra.MONARCHa.", cont, ".l250.2023-", mes, ".nc"))
  mean(rast(ruta))
}

# Aplicar a todos los meses y contaminantes
medias <- lapply(contaminantes, function(cont) {
  capas <- lapply(meses, media_mensual, cont = cont)
  mean(rast(capas))   # media anual
})

names(medias) <- contaminantes

# Acceder a los resultados
no2_anual  <- medias[["no2"]]
pm25_anual <- medias[["pm2p5"]]
# Verificar
print(no2_anual)
print(pm25_anual)

summary(no2_anual)# global(no2_anual,  "range")
summary(pm25_anual)# global(pm25_anual, "range")

# verificación visual
bbox_galicia <- st_bbox(municipios) # Bounding box para recortar rasters a Galcia
ext_gal <- ext(
  bbox_galicia["xmin"] - 0.5,
  bbox_galicia["xmax"] + 0.5,
  bbox_galicia["ymin"] - 0.5,
  bbox_galicia["ymax"] + 0.5
)

no2_mean  <- crop(no2_anual,  ext_gal)
pm25_mean <- crop(pm25_anual, ext_gal)
windows()
par(mfrow=c(2,1))
plot(no2_mean)
plot(pm25_mean)

#  3. Extracción Areal- Extraer el valor medio por municipio
municipios$no2_mean <- exact_extract(no2_mean, municipios, "mean")
municipios$pm25_mean <- exact_extract(pm25_mean, municipios, "mean")
head(st_drop_geometry(municipios[, c("name", "no2_mean", "pm25_mean")]))
summary(municipios$no2_mean) 
summary(municipios$pm25_mean)

# 4. Datos socioeconómicos: INE Atlas de Renta
# Renta primaria por hogar (GHI) según municipio
library(stringr)


GHI <- read.csv("datos/rentaprimariamun2023.csv", 
                encoding = "UTF-8")
ghi<-GHI[, 3:4]
colnames(ghi) <- c("municipio", "renta.hogar") # renta en miles de euros
head(ghi)

# 5. Limpieza de datos
# -- 5.1 Limpieza nombres de municipios en GHI (data.frame)
ghi<-ghi |> mutate(
  renta.hogar = as.numeric(renta.hogar),
  # Limpiar acentos y caracteres especiales 
  nombre_limpio = municipio |>
    str_to_lower() |>
    str_replace_all("á", "a") |>
    str_replace_all("é", "e") |>
    str_replace_all("í", "i") |>
    str_replace_all("ó", "o") |>
    str_replace_all("ú", "u") |>
    str_replace_all("ñ", "n") |>
    str_trim()
)
head(ghi)

# -- 5.2 Limpiar nombres en objeto municipios (sf)
municipios <- municipios |>
  mutate(
    nombre_limpio = name |>
      str_to_lower() |>
      str_replace_all("á", "a") |>
      str_replace_all("é", "e") |>
      str_replace_all("í", "i") |>
      str_replace_all("ó", "o") |>
      str_replace_all("ú", "u") |>
      str_replace_all("ñ", "n") |>
      str_trim()
  )

# -- 5.3 Verificar coincidencias de nombres
sin_match <- municipios |>
  st_drop_geometry() |>
  filter(!nombre_limpio %in% ghi$nombre_limpio) |>
  select(name, nombre_limpio)

cat("\nMunicipios sin match:\n")
print(sin_match)

# Corrección manual para que coincidan los nombres
ghi <- ghi |>
  mutate(nombre_limpio = case_when(
    nombre_limpio == "pereiro de aguiar"   ~ "pereiro de aguiar, o",
    TRUE ~ nombre_limpio           # el resto sin cambio
  ))


# Hacer el join por nombre limpio
municipios <- municipios |>
  left_join(
    ghi |> select(nombre_limpio, renta.hogar),
    by = "nombre_limpio"
  )

# Verificar si hay municipios sin datos
sum(is.na(municipios$renta.hogar))


# -- 5.2 Transformación logarítmica
municipios <- municipios |>
  mutate(
    renta_log = log(renta.hogar),
    no2_log = log(no2_mean),
    pm25_log = log(pm25_mean),
  )

summary(municipios$renta.hogar)
summary(municipios$renta_log)

# 6. Análisis Exploratorio espacial (ESDA)
# --6.1 Mapas coropléticos
# Visualizar extracción
windows()
library(tmap)
#tmap_mode("view")   # activa modo interactivo
# --- 5.1 Mapas coropléticos ---
tmap_mode("plot")

# NO2
mapa_no2 <- tm_shape(municipios) +
  tm_fill("no2_mean",
          style = "quantile", n = 5,
          palette = "YlOrRd",
          title = "NO2 medio (µg/m3)") +
  tm_borders(col = "white", lwd = 0.5) +
  tm_layout(title = "Concentración media anual de NO2 (2023)",
            legend.outside = TRUE)
#tmap_save(mapa_no2, "outputs/mapa_no2.png", dpi = 300)

# PM2.5
mapa_pm25 <- tm_shape(municipios) +
  tm_fill("pm25_mean",
          style = "quantile", n = 5,
          palette = "Purples",
          title = "PM2.5 medio (µg/m3)") +
  tm_borders(col = "white", lwd = 0.5) +
  tm_layout(title = "Concentración media anual de PM2.5 (2023)",
            legend.outside = TRUE)
#tmap_save(mapa_pm25, "outputs/mapa_pm25.png", dpi = 300)

# Renta media
mapa_renta <- tm_shape(municipios) +
  tm_fill("renta.hogar",
          style = "quantile", n = 5,
          palette = "Blues",
          title = "Renta primaria (euros/hogar)") +
  tm_borders(col = "white", lwd = 0.5) +
  tm_layout(title = "Renta neta media por hogar (2023)",
            legend.outside = TRUE)
#tmap_save(mapa_renta, "outputs/mapa_renta.png", dpi = 300)

# Panel comparativo
windows()
tmap_arrange(mapa_no2, mapa_pm25, mapa_renta, ncol = 3)

# -- 6.2 Estadísticas descriptivas
vars <- c("no2_mean", "pm25_mean", "renta.hogar")
summary(st_drop_geometry(provincias[, vars]))


# 7. Matriz de pesos o vecindario, W
