###
# Código R para solucionar el taller 1
###

# Librerías requeridas
library(ggplot2) # graficos
library(dplyr) # manipulación de datos
library(tidyr) # organizar los datos
library(xts) # manipulación de fechas

# Lectura de datos
datos_orig <- read.table(file = 'IndProductividad.txt', header = FALSE)
head(datos_orig)

# Organizamos los datos y le damos el formato adecuado
datos <- gather(as.data.frame(t(datos_orig)))[, 2]
datos <- xts(datos,
             as.yearmon(seq.Date(from = as.Date('1950-01-01', '%Y-%m-%d'),
                                 to = as.Date('1973-12-01', '%Y-%m-%d'),
                                 by = 'month')))
names(datos) <- 'ind_prod'
head(datos)

# Análisis exploratorio inicial
summary(datos)
plot(datos)

# Punto 1
decompose(datos)
