### ----------------------------------------------------------------------------
#
# Solución taller 4: técnicas de pronóstico
# Juan Pablo Calle Quintero
# jpcalleq@unal.edu.co
# 
### ----------------------------------------------------------------------------

## Librerías requeridas --------------------------------------------------------

library(data.table)
library(ggplot2)
library(gridExtra)
library(forecast)
library(car)
library(xtable)


## Definición funciones auxiliares ---------------------------------------------

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


## Lectura de datos ------------------------------------------------------------

datos_orig <- scan(file = '../datos/licor.txt', skip = 1)
datos_ts <- ts(datos_orig, start = c(1967, 1), end = c(1994, 12),
               frequency = 12)

## Organizamos los datos y le damos el formato adecuado
fecha <- seq.Date(from = as.Date('1967-01-01', '%Y-%m-%d'),
                  to = as.Date('1994-12-01', '%Y-%m-%d'),
                  by = 'month')

datos <- data.table(fecha, licor = datos_orig)
datos[, log_licor := log(licor)]
datos

## Análisis gráfico ------------------------------------------------------------

p1_base <- ggplot(data = datos, aes(x = fecha)) + theme_bw(10)

p1_serie <- p1_base + geom_line(aes(y = licor)) +
  labs(x = 'Tiempo',
       y = 'Ventas de licor [1000\'s USD]')
p1_log_serie <- p1_base + geom_line(aes(y = log_licor)) +
  labs(x = 'Tiempo',
       y = 'log(Ventas de licor [1000\'s USD])')

pdf('../img/p1_series.pdf', width = 7, height = 3.5)

grid.arrange(p1_serie, p1_log_serie, ncol = 2, nrow = 1)

dev.off()


## Punto 1 ---------------------------------------------------------------------

pdf('../img/p1_serie.pdf', width = 2.5, height = 2.3)

p1_base + geom_line(aes(y = licor), size = 0.1) + theme_bw(6) +
  labs(x = '', y = 'Ventas de licor [1000\'s USD]')

dev.off()


acf_serie <- acf(datos[, licor], lag.max = 48, plot = FALSE)$acf[-1]

pdf('../img/p1_acf.pdf', width = 2.5, height = 1.8)

opar <- par(no.readonly = TRUE)
par(mar = c(4, 4, 1, 1), cex.lab = 0.5, cex.axis = 0.5,
    tcl = -0.2, lwd = 0.6)

plot(acf_serie, ylim = c(-1, 1), type = 'h', xlab = 'k', ylab = 'ACF')
abline(h = 0)
text(12 * 1:4, 0.98, labels = 12 * 1:4, cex = 0.3)

par(opar)

dev.off()


## Punto 2 ---------------------------------------------------------------------

n_ajuste <- nrow(datos) - 12 # quitamos los últimos 12 datos para  ajustar
datos_ajuste <- datos[1:n_ajuste]
datos_ajuste[, t := 1:.N] # crear variable tiempo
ind_mes <- seasonaldummy(ts(datos_ajuste[, licor],
                            freq = 12, start = c(1967, 1)))
datos_ajuste[, I1 := ind_mes[, 1]]
datos_ajuste[, I2 := ind_mes[, 2]]
datos_ajuste[, I3 := ind_mes[, 3]]
datos_ajuste[, I4 := ind_mes[, 4]]
datos_ajuste[, I5 := ind_mes[, 5]]
datos_ajuste[, I6 := ind_mes[, 6]]
datos_ajuste[, I7 := ind_mes[, 7]]
datos_ajuste[, I8 := ind_mes[, 8]]
datos_ajuste[, I9 := ind_mes[, 9]]
datos_ajuste[, I10 := ind_mes[, 10]]
datos_ajuste[, I11 := ind_mes[, 11]]

mod1 <- lm(log_licor ~ t + I(t ^ 2) + I1 + I2 + I3 + I4 + I5 + I6 +
             I7 + I8 + I9 + I10 + I11, data = datos_ajuste)
summary(mod1)
xtable(summary(mod1),
       caption = 'Parámetros estimados modelo 1 log-cuadrático estacional.',
       label = 'tab:p2_mod1_ajuste',
       align = 'lrrrr')

val_ini <- coef(mod1) # valores iniciales modelo exponencial cuadrático

mod2 <- nls(licor ~ exp(b0 + b1*t + b2*I(t^2) +
                       d1*I1 + d2*I2 + d3*I3 + d4*I4 + d5*I5 + d6*I6 + d7*I7 +
                       d8*I8 + d9*I9 + d10*I10 + d11*I11),
            start = list(b0 = val_ini[1], b1 = val_ini[2], b2 = val_ini[3],
                         d1 = val_ini[4], d2 = val_ini[5], d3 = val_ini[6],
                         d4 = val_ini[7], d5 = val_ini[8], d6 = val_ini[9],
                         d7 = val_ini[10], d8 = val_ini[11], d9 = val_ini[12],
                         d10 = val_ini[13], d11 = val_ini[14]),
            data = datos_ajuste)
summary(mod2)

datos_ajuste[, res_mod1 := residuals(mod1)]
datos_ajuste[, res_mod2 := residuals(mod2)]
datos_ajuste[, aju_mod1 := fitted(mod1)]
datos_ajuste[, aju_mod2 := fitted(mod2)]

# Gráficos de residuales
p23_base <- ggplot(data = datos_ajuste) + theme_bw(8)

p2_mod1_tpo <- p23_base + geom_line(aes(x = fecha, y = res_mod1),
                                    size = 0.5) +
  labs(x = 'Tiempo', y = 'Residuales')

p2_mod1_aju <- p23_base + geom_point(aes(x = aju_mod1, y = res_mod1),
                                    size = 2) +
  labs(x = 'Ajustados', y = 'Residuales')

p2_mod2_tpo <- p23_base + geom_line(aes(x = fecha, y = res_mod2),
                                    size = 0.5) +
  labs(x = 'Tiempo', y = 'Residuales')

p2_mod2_aju <- p23_base + geom_point(aes(x = aju_mod2, y = res_mod2),
                                     size = 2) +
  labs(x = 'Ajustados', y = 'Residuales')

pdf('../img/p2_residuales_mod1.pdf', height = 2.5, width = 7)

grid.arrange(p2_mod1_tpo, p2_mod1_aju, nrow = 1, ncol = 2)

dev.off()

pdf('../img/p2_residuales_mod2.pdf', height = 2.5, width = 7)

grid.arrange(p2_mod2_tpo, p2_mod2_aju, nrow = 1, ncol = 2)

dev.off()

# Correlogramas muestrales de los residuales

acf_res_mod1 <- acf(datos_ajuste$res_mod1, lag.max = 48, plot = FALSE)$acf[-1]
acf_res_mod2 <- acf(datos_ajuste$res_mod2, lag.max = 48, plot = FALSE)$acf[-1]
pacf_res_mod1 <- pacf(datos_ajuste$res_mod1, lag.max = 48, plot = FALSE)$acf
pacf_res_mod2 <- pacf(datos_ajuste$res_mod2, lag.max = 48, plot = FALSE)$acf


pdf('../img/p2_correlogramas_mod1.pdf', height = 2, width = 7)

opar <- par(no.readonly = TRUE)
par(mar = c(4, 4, 1, 1), cex.lab = 0.5, cex.axis = 0.5,
    tcl = -0.2, lwd = 0.6)
layout(matrix(c(1, 2), 1, 2))

acf(datos_ajuste$res_mod1, lag.max = 48, ci.type = 'ma',
    ylim = c(-1, 1), xlab = 'Rezago', ylab = 'ACF')
pacf(datos_ajuste$res_mod1, lag.max = 48,
     ylim = c(-1, 1), xlab = 'Rezago', ylab = 'ACF')

layout(matrix(c(1, 1), 1, 1))
par(opar)

dev.off()


pdf('../img/p2_correlogramas_mod2.pdf', height = 2, width = 7)

opar <- par(no.readonly = TRUE)
par(mar = c(4, 4, 1, 1), cex.lab = 0.5, cex.axis = 0.5,
    tcl = -0.2, lwd = 0.6)
layout(matrix(c(1, 2), 1, 2))

acf(datos_ajuste$res_mod2, lag.max = 48, ci.type = 'ma',
    ylim = c(-1, 1), xlab = 'Rezago', ylab = 'ACF')
pacf(datos_ajuste$res_mod2, lag.max = 48,
     ylim = c(-1, 1), xlab = 'Rezago', ylab = 'ACF')

layout(matrix(c(1, 1), 1, 1))
par(opar)

dev.off()

xtable(
pruebaBPJB(datos_ajuste[, res_mod1], maxlag = 48, type = 'Box'),
caption = 'Resultados prueba Box-Pierce modelo 1.',
label = 'tab:p2_bp'
)
xtable(
  pruebaBPJB(datos_ajuste[, res_mod1], maxlag = 48, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box modelo 1.',
  label = 'tab:p2_lb'
)

xtable(
pruebaBPJB(datos_ajuste[, res_mod2], maxlag = 48, type = 'Box'),
caption = 'Resultados prueba Box-Pierce modelo 2.',
label = 'tab:p3_bp'
)

xtable(
pruebaBPJB(datos_ajuste[, res_mod2], maxlag = 48, type = 'Ljung'),
caption = 'Resultados prueba Ljung-Box modelo 2.',
label = 'tab:p3_lb'
)

xtable(
pruebaDW(mod1),
caption = 'Resultados prueba de orden 1 Durbin-Watson.',
label = 'tab:p2_db'
)

pruebaDW(mod2) # error, el modelo no es lineal
