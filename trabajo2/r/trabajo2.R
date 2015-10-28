### ----------------------------------------------------------------------------
#
# Trabajo 2: Ajuste de modelos de regresión con errores estructurales ARMA
# Juan Pablo Calle Quintero
# jpcalleq@unal.edu.co
# 
### ----------------------------------------------------------------------------


## Librerías requeridas --------------------------------------------------------

library(ggplot2) # gráficos de calidad
library(grid) # para poder utilizar la función unit() en ggplot2
library(gridExtra) # utilidades adicionales para graficar grid.arrange()
library(scales) # formatear los valores de los ejes de las gráficas.
library(data.table) # manipulación de datos
library(xtable) # exportar objetos a LaTeX
library(forecast)
library(car)


## Definición de funciones auxiliares ------------------------------------------

plot_pacf <- function(x, lag = NULL, ciACF = FALSE, ciPACF = FALSE){
  # Graficar ACF y PACF de una serie de tiempo
  
  op <- par(no.readonly = TRUE)
  
  if (is.null(lag)) {
    lag <- round(12 * (length(x) / 100) ^ (1 / 4), 0)
  }
  
  U <- 2 / sqrt(length(x))
  L <- -U
  
  dacf <- acf(x = x, lag.max = lag, plot = FALSE)$acf[-1]
  
  layout(matrix(c(1, 2), 2, 1))
  par(mar = c(0.5, 4, 4, 2) + 0.1, lab = c(lag, 5 , 7))
  plot(x = c(1:lag), y = dacf, ylim = c(-1, 1), xlab = "", ylab = "ACF",
       type = "h", axes = FALSE, main = "Correlogramas", xlim = c(1, lag))
  box(col = "gray")
  abline(h = 0, lty = 1)
  if (ciACF) {
    abline(h = c(L, U), lty = c(2, 2), col = c(4, 4))
  }
  axis(side = 2, col = "gray", cex.axis = 0.5)
  
  dpacf <- pacf(x = x, lag.max = lag, plot = FALSE)$acf
  
  par(mar = c(4, 4, 0.5, 2) + 0.1, lab = c(lag, 5 , 7))
  plot(x = c(1:lag), y = dpacf, ylim = c(-1, 1), xlab = "Lag", ylab = "PACF",
       type = "h", axes = FALSE, xlim = c(1, lag))
  box(col = "gray")
  abline(h = 0, lty = 1)
  if (ciPACF) {
    abline(h = c(L, U), lty = c(2, 2), col = c(4, 4))
  }
  axis(side = 2, col = "gray", cex.axis = 0.5)
  axis(side = 1, col = "gray", cex.axis = 0.5)
  
  par(op)
}

# Gráfico en blanco, a veces útil para usar con la función grid.arrange()
blankPlot <- ggplot() + geom_blank(aes(1, 1)) +
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pruebaDW <- function(modelo) {
  # Calcular test de Durbin-Watson
  # modelos: un objeto de tipo modelo lineal
  
  dwneg <- durbinWatsonTest(modelo, max.lag = 1, method = "normal",
                            alternative = "negative")
  
  dwpos <- durbinWatsonTest(modelo, max.lag = 1, method = "normal",
                            alternative = "positive")
  
  res <- data.frame(1, dwneg$r, dwneg$dw, dwpos$p, dwneg$p)
  
  names(res) <- c("lag", "rho estimado", "Estadístico D-W", "VP rho>0",
                  "VP rho<0")
  
  return(res)
}

pruebaBPJB <- function(serie, maxlag, type = "Box") {
  # Calcular test Box-Pierce y Ljung-Box
  
  aux <- floor(maxlag / 6)
  X.squared <- c(rep(NA, aux))
  df <- c(rep(NA, aux))
  p.value <- c(rep(NA, aux))
  
  for (i in seq_len(aux)) {
    test <- Box.test(serie, lag = 6 * i, type = type)
    X.squared[i] <- test[[1]]
    df[i] <- test[[2]]
    p.value[i] <- test[[3]]
  }
  
  lag <- 6 * c(1:aux)
  test_res <- as.data.frame(cbind(X.squared, df, p.value))
  rownames(test_res) <- lag
  
  return(test_res)
}

C_p <- function(residuales, type = 'AIC', p = 1){
  # calcula criterio de comparación de modelos C*(p)
  # recibe los residuales en un vector y el tipo de criterio, AIC o BIA
  # p = número de parámetros, incluyendo el intercepto
  
  n <- length(residuales)
  
  if (type == 'AIC'){
    c_p <- log(mean(residuales^2)) + 2 * p/n
  } else if (type == 'BIC'){
    c_p <- log(mean(residuales^2)) + log(n) * p/n
  }
  
  return(c_p)
}


## Lectura de datos ------------------------------------------------------------

datos_orig <- fread('../datos/datos.csv', header = FALSE, skip = 11,
                         sep = ';', dec = ',',
                         colClasses = c(rep("NULL", 12),
                                        "numeric",
                                        rep("NULL", 5)))
setnames(datos_orig, 'ivand')
datos_orig[, log_ivand := log(ivand)]

nro_datos <- datos_orig[, .N]
datos_orig_ts <- ts(datos_orig$ivand, freq = 4, start = c(2000, 1))
datos_log_ts <- ts(datos_orig$log_ivand, freq = 4, start = c(2000, 1))

# organizamos los datos en undata frame y creamos otras variables de interés
datos_orig[,
           ':='(fecha = seq.Date(from = as.Date('2000-03-01', '%Y-%m-%d'),
                                 length.out = nro_datos,
                                 by = 'quarter'),
                tri = gl(n = 4, k = 1, length = nro_datos,
                         labels = c('T1', 'T2', 'T3', 'T4')),
                t = 1:nro_datos
           )
           ]


## Análisis exploratorio inicial -----------------------------------------------
summary(datos_orig)
plot(datos_orig_ts)
abline(v = 2000:2015, lty = 2)
plot(datos_log_ts)
abline(v = 2000:2015, lty = 2)

## Punto 2 ---------------------------------------------------------------------

# Gráficos de las series; original y su logaritmo
serie_orig <- ggplot(data = datos_orig, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand)) +
  labs(x = '',
       y = 'IVA no deducible [1000\'s mill pesos]') +
  scale_y_continuous(labels = dollar)
serie_orig

serie_log <- ggplot(data = datos_orig, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = log_ivand)) +
  labs(x = '',
       y = 'log(IVA no deducible [1000\'s mill pesos])')
serie_log

pdf(file = '../img/series.pdf', width = 7, height = 3.5)

grid.arrange(serie_orig, serie_log, nrow = 1, ncol = 2)

dev.off()

# ACF del logaritmo de la serie
pdf(file = '../img/acf_log_serie.pdf', width = 6, height = 4)

datos_ajuste[, acf(log_ivand, ci.type = 'ma',
                   ylim = c(-1, 1), xlim = c(1, 10 * log10(.N)),
                   main = '', xlab = 'k')]

dev.off()

## Ajuste de modelos globales --------------------------------------------------

datos_ajuste <- datos_orig[1:(nro_datos - 4)]
datos_ajuste[, c('I2', 'I3', 'I4') := list(
  ifelse(tri == 'T2', 1, 0),
  ifelse(tri == 'T3', 1, 0),
  ifelse(tri == 'T4', 1, 0)
)]

# modelo auxiliar para estimar modelo no lineal
mod_aux <- lm(log_ivand ~ t + I(t ^ 2) + I2 + I3 + I4, data = datos_ajuste)
summary(mod_aux)

# modelo: exponencial cuadrático estacional
val_ini <- coef(mod_aux)

modelo <- nls(ivand ~ exp(b0 + b1 * t + b2 * I(t ^ 2) +
                            d2 * I2 + d3 * I3 + d4 * I4),
            start = list(b0 = val_ini[1], b1 = val_ini[2], b2 = val_ini[3],
                         d2 = val_ini[4], d3 = val_ini[5], d4 = val_ini[6]),
            data = datos_ajuste)

summary(modelo)

datos_ajuste[, aju := fitted(modelo)] # valores ajustados
datos_ajuste[, res := residuals(modelo)] # error E_t


## Postulación de modeloa ARMA para el error E_t -------------------------------

# Gráficos de los residuales
res_vs_t <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res)) +
  geom_hline(yintercept = 0, colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = 'Residuales')

res_vs_aju <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju, y = res)) +
  geom_hline(yintercept = 0, colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = 'Residuales')

pdf(file = '../img/residuales.pdf', width = 7, height = 3.5)

grid.arrange(res_vs_t, res_vs_aju, nrow = 1, ncol = 2)

dev.off()

# ACF y PACF de los errores
pdf(file = '../img/acf_pacf_res.pdf', width = 4, height = 6)

opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))

datos_ajuste[, acf(res, ci.type = 'w',
                   ylim = c(-1, 1), xlim = c(1, 10 * log10(.N)),
                   main = '', xlab = 'k')]

datos_ajuste[, pacf(res,
                   ylim = c(-1, 1), xlim = c(1, 10 * log10(.N)),
                   main = '', xlab = 'k', ylab = 'PACF')]

par(opar)

dev.off()
