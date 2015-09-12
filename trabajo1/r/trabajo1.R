#------------------------------------------------------------------------------#
# Código R trabajo 1: Análisis PIB Colombia por Ramas Actividad Económica      #
# Autor: Juan Pablo Calle Quintero                                             #
# Fecha: 24/09/2015                                                            #
# Correo: jpcalleq@unal.edu.co                                                 #
#------------------------------------------------------------------------------#

## Librerías requeridas --------------------------------------------------------
library(ggplot2) # graficos de calidad
library(grid) # para poder utilizar la función unit() en ggplot2
library(gridExtra) # funciones auxiliares graficación múltiple, grid.arrange()
library(dplyr) # manipulación de datos
library(xtable) # exportar tablas en formato LaTeX
library(scales) # formatear los valores de los ejes de las gráficas.
library(forecast) # funciones de pronóstico para series de tiempo y mod lineales
# theme_set(theme_bw(base_size = 7)) # definir tamaño de fuente gráficos ggplot2

## Funciones auxiliares --------------------------------------------------------
ggqqplot <- function(mod, resid.type = 'student'){
  # Devuelve un gráfico qqplot al estilo ggplot2
  
  if (resid.type == 'student'){
    res_std <- rstudent(mod)
  } else if (resid.type == 'regular'){
    # cuando el modelo ajustado es o lineal (función nls())
    res_std <- resid(mod) / sd(resid(mod))
  }
  
  val_aju <- fitted(mod)
  
  # Almacenar los datos en un data frame para poder usar ggplot
  dat_tmp <- data.frame(res_std, val_aju)
  
  qq_val <- qqnorm(res_std, plot.it = FALSE)
  dat_tmp$q_teo <- qq_val$x # añadir cuantiles teóricos al data frame
  dat_tmp$q_muestra <- qq_val$y # añadir cuantiles muestrales al data frame
  
  slope_qq <- diff(quantile(res_std, c(0.25, 0.75), type = 7)) /
    diff(qnorm(c(0.25, 0.75)))
  
  inter_qq <- quantile(res_std, c(0.25, 0.75), type = 7)[1L] -
    slope_qq * qnorm(c(0.25, 0.75))[1L]
  
  # pruebas de normalidad Shapiro-Wilk y Kolmogorov-Smirnov
  sw_vp <- shapiro.test(res_std)$p.value
  ks_vp <- ks.test(res_std, 'pnorm')$p.value
  
  pruebas_norm <- paste('Prueba  p-valor\n',
                        'SW        ', round(sw_vp, 4), '\n',
                        'KS         ', round(ks_vp, 4), sep = '')
  
  ub_y <- max(dat_tmp$q_muestra) - 0.2 # ubicación de los p-valores
  ub_x <- min(dat_tmp$q_teo)
  
  # qqplot con ggplot
  ggplot(data = dat_tmp, aes(x = q_teo, y = q_muestra)) + theme_bw(8) +
    geom_point(size = 1) +
    geom_abline(intercept = inter_qq, y = slope_qq, colour = 'red') +
    labs(x = 'Cuantiles Teóricos', y = 'Cuantiles Muestrales',
         title = 'QQ-Plot') +
    annotate('text', x = ub_x, y = ub_y, size = 2, parse = FALSE,
             label = pruebas_norm, hjust = 0) +
    theme(plot.margin = unit(c(0.1,0.1,0,0), "lines"),
          panel.border = element_blank(),
          axis.ticks = element_blank())
}

C_p <- function(residuales, type = 'AIC', p = 1){
  # calcula criterio de comparación de modelos C*(p)
  # recibe los residuales en un vector y el tipo de criterio, AIC o BIA
  # p = número de parámetros incluido el intercepto
  
  n <- length(residuales)
  
  if (type == 'AIC'){
    c_p <- log(mean(residuales^2)) + 2 * p/n
  } else if (type == 'BIC'){
    c_p <- log(mean(residuales^2)) + log(n) * p/n
  }
  
  return(c_p)
}


## Lectura de datos ------------------------------------------------------------
datos_orig <- read.table('../datos/datos.csv', header = FALSE, skip = 11,
                         sep=';', dec = ',',
                         colClasses = c(rep("NULL", 12),
                                        "numeric",
                                        rep("NULL", 5)))
names(datos_orig) <- 'pib'
datos_orig_ts <- ts(datos_orig, freq = 4, start = c(2000, 1))
datos_log_ts <- ts(log(datos_orig), freq = 4, start = c(2000, 1))
datos_orig
datos_orig_ts
datos_log_ts

## Organizamos los datos y le damos el formato adecuado
fecha <- seq.Date(from = as.Date('2000-03-01', '%Y-%m-%d'),
                  length.out = 59,
                  by = 'quarter')
mes <- reorder(months.Date(fecha), c(rep(c(3, 6, 9, 12), 14), 3, 6, 9))
relevel(mes, ref = 'December')
log_pib <- log(datos_orig)
names(log_pib) <- 'log_pib'
datos <- data.frame(mes, fecha, datos_orig, log_pib)
head(datos)
str(datos)

## Análisis exploratorio inicial -----------------------------------------------
summary(datos)
plot(datos_orig_ts)


## 1. Análisis descriptivo de la serie -----------------------------------------

# Gráficos de la serie
p1_serie_orig <- ggplot(data = datos, aes(x = fecha)) +
  theme_bw() +
  geom_line(aes(y = pib)) +
  labs(x = '',
       y = 'PIB [1000\'s mill pesos]')
p1_serie_orig

p1_serie_log <- ggplot(data = datos, aes(x = fecha)) +
  theme_bw() +
  geom_line(aes(y = log_pib)) +
  labs(x = '',
       y = 'log(PIB [1000\'s mill pesos])')
p1_serie_log


## 2. Postulación de modelos ---------------------------------------------------

# Gráficos descomposición de la serie
ts_desc_orig <- decompose(datos_orig_ts)

pdf('../img/p1_descomp.pdf', width = 7, height = 4)
opar <- par()
par(mfrow = c(4, 1),
    oma = c(1, 1, 1, 1),
    mar = c(1, 4, 1, 1),
    xaxt = 'n',
    cex.axis = 1,
    cex.lab = 1,
    lwd = 0.5)
plot(ts_desc_orig$x,
     xlab = '',
     ylab = expression(paste(Y[t])))
plot(ts_desc_orig$trend,
     xlab = '',
     ylab = expression(paste(T[t])))
plot(ts_desc_orig$seasonal,
     xlab = '',
     ylab = expression(paste(S[t])))
abline(v = 2000:2015, lty = 2, col =  'red')
par(xaxt = 's')
plot(ts_desc_orig$random,
     xlab = 'Tiempo',
     ylab = expression(paste(E[t])))
abline(h = 0, lty = 2, col =  'red')
par(opar)
dev.off()

ts_desc_log <- decompose(datos_log_ts, type = 'additive')

pdf('../img/p1_descomp_log.pdf', width = 7, height = 4)
opar <- par()
par(mfrow = c(4, 1),
    oma = c(1, 1, 1, 1),
    mar = c(1, 4, 1, 1),
    xaxt = 'n',
    cex.axis = 1,
    cex.lab = 1,
    lwd = 0.5)
plot(ts_desc_log$x,
     xlab = '',
     ylab = expression(paste(Y[t])))
plot(ts_desc_log$trend,
     xlab = '',
     ylab = expression(paste(T[t])))
plot(ts_desc_log$seasonal,
     xlab = '',
     ylab = expression(paste(S[t])))
abline(v = 2000:2015, lty = 2, col =  'red')
par(xaxt = 's')
plot(ts_desc_log$random,
     xlab = 'Tiempo',
     ylab = expression(paste(E[t])))
abline(h = 0, lty = 2, col =  'red')
par(opar)
dev.off()

# Box plot
prom_trim_orig <- summarise(group_by(datos, mes), prom = mean(pib))
prom_trim_log <- summarise(group_by(datos, mes), prom = mean(log_pib))

trimestres <- c('T1', 'T2', 'T3', 'T4')


p1_boxplot_mes_orig <- ggplot(data = datos) +
  theme_bw() +
  geom_boxplot(aes(x = mes, y = pib)) +
  geom_line(data = prom_trim_orig, aes(x = 1:4, y = prom),
            size = 0.2, colour = 'red', linetype = 6) +
  scale_x_discrete(labels = trimestres) +
  labs(x = '', y = 'PIB [1000\'s mill pesos]') +
  theme(
    axis.line = element_blank(),
    #panel.border = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_y_continuous(labels = dollar)

p1_boxplot_mes_orig

p1_boxplot_mes_log <- ggplot(data = datos) +
  theme_bw() +
  geom_boxplot(aes(x = mes, y = log_pib)) +
  geom_line(data = prom_mes_log, aes(x = 1:4, y = prom),
            size = 0.2, colour = 'red', linetype = 6) +
  scale_x_discrete(labels = trimestres) +
  labs(x = '', y = 'log(PIB [1000\'s mill pesos])') +
  theme(
    axis.line = element_blank(),
    #panel.border = element_blank(),
    axis.ticks = element_blank()
  )

p1_boxplot_mes_log

# Periodogramas
diff_pib <- diff(datos$pib)
diff_log_pib <- diff(datos$log_pib)

nro_datos_dif <- length(diff_pib)

period_pib <- abs(fft(diff_pib) / sqrt(nro_datos_dif))^2
period_log_pib <- abs(fft(diff_log_pib) / sqrt(nro_datos_dif))^2

period_pib_scale <- period_pib * (4/nro_datos_dif)
period_log_pib_scale <- period_log_pib * (4/nro_datos_dif)

freq <- (0:(nro_datos_dif - 1)) / nro_datos_dif
plot(freq, period_pib_scale, type = 'l')
plot(freq, period_log_pib_scale, type = 'l')

periodogramas <- data.frame(freq,
                            pib = period_pib_scale,
                            log_pib = period_log_pib_scale)

p1_periodograma_pib <- ggplot(data = periodogramas) + theme_bw() +
  geom_line(aes(x = freq, y = pib)) +
  geom_vline(xintercept = (1:3)/4, linetype = 2, colour = 'red') +
  labs(y = 'Periodograma', x = 'Frecuencia') +
  scale_x_continuous(limits = c(0, 0.5))

p1_periodograma_pib

p1_periodograma_log_pib <- ggplot(data = periodogramas) + theme_bw() +
  geom_line(aes(x = freq, y = log_pib)) +
  geom_vline(xintercept = (1:3)/4, linetype = 2, colour = 'red') +
  labs(y = 'Periodograma', x = 'Frecuencia') +
  scale_x_continuous(limits = c(0, 0.5))

p1_periodograma_log_pib


## 3. Ajuste de modelos con validación cruzada ---------------------------------


## 4. Análisis de residuales y validación de supuestos -------------------------


## 5. Pronósticos para validación cruzada --------------------------------------


## 6. Selección del mejor modelo -----------------------------------------------
