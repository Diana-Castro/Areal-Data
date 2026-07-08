# TFM: ANÁLISIS DE DATOS AGREGADOS POR ÁREA- contaminación atmosférica y desigualdad socioeconómica
# Preprocesado de datos y creación de base de datos final para análisis espacial
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

# 1. Límites Administrativos: municipios de Galicia
municipios <- st_transform(esp_get_munic(region = "ES11"), crs= 4326) 
# municipios <- esp_get_munic() |> filter(codauto == "12") |> st_transform(crs = 4326)
# Verificar
print(municipios)
windows()
plot(st_geometry(municipios), main = "Municipios de Galicia")
# mapa con leaflet
leaflet(municipios) |> 
  addTiles() |> 
  addPolygons(color = "blue", weight = 1, fillOpacity = 0.5)
 
# 2. Cargar el raster NetCDF de datos CAMS
meses <- sprintf("%02d", 1:12) 
contaminantes <- c("no2", "pm2p5")

# Función que calcula la mediana mensual
mediana_mensual <- function(mes, cont) {
  ruta <- file.path("data",paste0("cams.eaq.vra.MONARCHa.", cont, ".l100.2023-", mes, ".nc"))
  median(rast(ruta)) #
}

# Aplicar a todos los meses y contaminantes
medianas <- lapply(contaminantes, function(cont) {
  capas <- lapply(meses, mediana_mensual, cont = cont)
  median(rast(capas))   # mediana anual
})

names(medianas) <- contaminantes

# Acceder a los resultados
no2_anual  <- medianas[["no2"]]
pm25_anual <- medianas[["pm2p5"]]
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

no2_median  <- crop(no2_anual,  ext_gal)
pm25_median <- crop(pm25_anual, ext_gal)
windows()
par(mfrow=c(2,1))
plot(no2_median)
plot(pm25_median)

#  3. Extracción Areal- Extraer el valor de mediana por municipio
municipios$no2_median <- exact_extract(no2_median, municipios, "median")
municipios$pm25_median <- exact_extract(pm25_median, municipios, "median")
head(st_drop_geometry(municipios[, c("name", "no2_median", "pm25_median")]))
summary(municipios$no2_median) 
summary(municipios$pm25_median)

# 4. Datos socioeconómicos: INE Atlas de Renta
# Renta per capita (pci) según municipio
library(stringr)


PCI <- read.csv("data/rentapcapmun2023.csv", 
                encoding = "UTF-8")
pci <-PCI[, 6:7]
colnames(pci) <- c("municipio", "renta.pc") # renta en miles de euros
head(pci)

# 5. Limpieza de datos
# -- 5.1 Limpieza nombres de municipios en PCI (data.frame)
pci<-pci |> mutate(
  renta.pc = as.numeric(renta.pc),
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
head(pci)

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
  filter(!nombre_limpio %in% pci$nombre_limpio) |>
  select(name, nombre_limpio)

cat("\nMunicipios sin match:\n")
print(sin_match)

# Corrección manual para que coincidan los nombres
pci <- pci |>
  mutate(nombre_limpio = case_when(
    nombre_limpio == "pereiro de aguiar"   ~ "pereiro de aguiar, o",
    TRUE ~ nombre_limpio           # el resto sin cambio
  ))


# Hacer el join por nombre limpio
municipios <- municipios |>
  left_join(
    pci |> select(nombre_limpio, renta.pc),
    by = "nombre_limpio"
  )

# Verificar si hay municipios sin datos
sum(is.na(municipios$renta.pc))


# -- 5.4 Transformación logarítmica
municipios <- municipios |>
  mutate(
    renta_log = log(renta.pc),
    no2_log = log(no2_median),
    pm25_log = log(pm25_median),
  )

summary(municipios$renta.pc)
summary(municipios$renta_log)

# Visualización de extracción

#tmap_mode("view")   # activa modo interactivo
# --- 5.1 Mapas coropléticos ---
tmap_mode("plot")

# NO2
mapa_no2 <- tm_shape(municipios) +
  tm_polygons("no2_median",
              style = "quantile",# "fixed"
              n=5,#breaks = c(-Inf, 10, 20, 30, 40, Inf),
              #labels = c("< 10", "10-20", "21–30", "31–40", "> 40"),
              palette = "brewer.yl_or_rd",
              title = "NO2 medio (µg/m3)",
              border.col = "white",
              border.alpha = 0.5) +
  tm_layout(title = "Concentración media anual de NO2 (2023)",
            legend.outside = TRUE)

# PM2.5
mapa_pm25 <- tm_shape(municipios) +
  tm_polygons("pm25_median",
              style = "quantile",# "fixed"
              n=5,#breaks = c(-Inf, 5, 10, 15, 20, 25, Inf),
              #labels = c("< 5", "5–10", "11–15", "16–20", "21–25", "> 25"),
              palette = "brewer.purples",
              title = "PM2.5 medio (µg/m3)",
              border.col = "white",
              border.alpha = 0.5) +
  tm_layout(title = "Concentración media anual de PM2.5 (2023)",
            legend.outside = TRUE)

# Renta media
mapa_renta <- tm_shape(municipios) +
  tm_fill("renta.pc",
          style = "quantile", n = 5,
          palette = "brewer.blues",
          title = "Renta bruta (miles de euros per cápita)") +
  tm_borders(col = "white", lwd = 0.5) +
  tm_layout(title = "Renta primaria bruta per cápita (2023)",
            legend.outside = TRUE)

# Panel comparativo
windows()
tmap_arrange(mapa_no2, mapa_pm25, mapa_renta, ncol = 3)


###### Datos limpios
# 5.5 Guardar objeto sf con datos limpios
st_write(municipios, "data/galmun_contaminacion_renta.gpkg", delete_dsn = FALSE) 
st_write(municipios, "data/galmun_contaminacion_renta.csv", delete_dsn = FALSE)


