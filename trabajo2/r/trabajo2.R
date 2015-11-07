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
library(forecast) # ajustar modelos ARIMA Arima()
library(TSA) # función eacf() y armasubsets()
library(car) # funcione DurbinWatsonTest
library(FitAR) # función SelectModel()


## Definición de funciones auxiliares ------------------------------------------

plot_pacf <- function(x, lag = NULL, ciACF = FALSE, ciPACF = FALSE){
  # Graficar ACF y PACF de una serie de tiempo
  
  op <- par(no.readonly = TRUE)
  
  if (is.null(lag)) {
    lag <- round(12 * (length(x) / 100) ^ (1 / 4), 0)
  }
  
  U <- 2 / sqrt(length(x))
  L <- -U
  
  dacf <- stats::acf(x = x, lag.max = lag, plot = FALSE)$acf[-1]
  
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
  # modelo: un objeto de tipo modelo lineal
  
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
  
  if (type == 'AIC') {
    c_p <- log(mean(residuales ^ 2)) + 2 * p/n
  } else if (type == 'BIC') {
    c_p <- log(mean(residuales ^ 2)) + log(n) * p/n
  }
  
  return(c_p)
}

ggqqplot <- function(mod) {
  # qqplot al estilo ggplot2
  
  residuales <- resid(mod)
  
  mu_resid <- mean(residuales)
  sd_resid <- sd(residuales)
  
  residuales <- resid(mod) / sd_resid
  
  # Almacenar los datos en un data frame para poder usar ggplot
  dat_tmp <- data.frame(residuales)
  
  qq_val <- qqnorm(residuales, plot.it = FALSE)
  dat_tmp$q_teo <- qq_val$x # añadir cuantiles teóricos al data frame
  dat_tmp$q_muestra <- qq_val$y # añadir cuantiles muestrales al data frame
  
  slope_qq <- diff(quantile(residuales, c(0.25, 0.75), type = 7)) /
    diff(qnorm(c(0.25, 0.75), mean = mu_resid, sd = sd_resid))
  
  inter_qq <- quantile(residuales, c(0.25, 0.75), type = 7)[1L] -
    slope_qq * qnorm(c(0.25, 0.75), mean = mu_resid, sd = sd_resid)[1L]
  
  # pruebas de normalidad Shapiro-Wilk y Kolmogorov-Smirnov
  sw_vp <- prettyNum(shapiro.test(residuales)$p.value)
  ks_vp <- prettyNum(ks.test(residuales, 'pnorm')$p.value)
  
  pruebas_norm <- paste('    valor p\n',
                        'SW: ', sw_vp, '\n',
                        'KS: ', ks_vp, sep = '')
  
  ub_y <- max(dat_tmp$q_muestra) - 0.5 # ubicación de los p-valores
  ub_x <- min(dat_tmp$q_teo) + 0.3
  
  # qqplot con ggplot
  ggplot(data = dat_tmp, aes(x = q_teo, y = q_muestra)) + theme_bw(10) +
    geom_point(size = 1) +
    geom_abline(intercept = inter_qq, y = slope_qq, colour = 'red') +
    labs(x = 'Cuantiles Teóricos', y = 'Cuantiles Muestrales') +
    annotate('text', x = ub_x, y = ub_y, size = 3, parse = FALSE,
             label = pruebas_norm, hjust = 0, family = 'Courier')
}

resumen_arima <- function(mod) {
  # mod: un modelos Arima()
  # Devuelve un resumen de la estimación del modelo
  
  df_mod <- mod$nobs - length(mod$coef[mod$mask])
  t0 <- mod$coef[mod$mask] / sqrt(diag(mod$var.coef))
  val_p <- pt(abs(t0), df = df_mod, lower.tail = FALSE) +
    pt(-abs(t0), df = df_mod)
  
  res <- cbind(mod$coef[mod$mask], sqrt(diag(mod$var.coef)), val_p)
  
  colnames(res) <- c('Estimación',
                     'Error Estándar',
                     paste('$Pr(|t_{', df_mod, '}| > |t_0|$', sep = ''))
  
  return(res)
}

calc_coef_recursivos <- function(modelo, k.ini = length(coef(modelo)) + 1) {
  # calcula los coeficientes del modelo recursivamente
  
  if (is(modelo, 'lm')) {
    
    npar_mod <- length(coef(modelo))
    nro_datos <- nrow(modelo$model)
    nro_obs_mod_i <- k.ini:nro_datos
    
    # La matriz para almacenar los valores de los coeficientes tiene en la
    # primera columna el número de observaciones usadas. En las siguientes
    # se almacenan los valores de los coeficientes con sus interv de conf
    mtrx_coef <- matrix(data = NA,
                        nrow = length(nro_obs_mod_i),
                        ncol = npar_mod * 3 + 1)
    
    idx_coef <- seq(2, ncol(mtrx_coef), 3)
    idx_lic <- seq(3, ncol(mtrx_coef), 3)
    idx_lsc <- seq(4, ncol(mtrx_coef), 3)
    
    colnames(mtrx_coef) <- paste('col', 1:ncol(mtrx_coef), sep = '.')
    colnames(mtrx_coef)[1] <- 'k'
    colnames(mtrx_coef)[idx_coef] <- names(coef(modelo))
    colnames(mtrx_coef)[idx_lic] <- paste(names(coef(modelo)), 'lic', sep = '.')
    colnames(mtrx_coef)[idx_lsc] <- paste(names(coef(modelo)), 'lsc', sep = '.')
    
    for (i in nro_obs_mod_i) {
      mod_tmp <- lm(formula(modelo),
                    data = modelo$model[1:i, ])
      sum_tmp <- summary(mod_tmp)
      coef_tmp <- coef(mod_tmp)
      t_0 <- qt(0.025, mod_tmp$df.residual)
      lic <- coef_tmp - t_0 * sum_tmp$coefficients[, 2]
      lsc <- coef_tmp + t_0 * sum_tmp$coefficients[, 2]
      mtrx_coef[i - k_ini + 1, 1] <- i
      mtrx_coef[i - k_ini + 1, idx_coef] <- coef_tmp
      mtrx_coef[i - k_ini + 1, idx_lic] <- lic
      mtrx_coef[i - k_ini + 1, idx_lsc] <- lsc
    }
    
  } else if (is(modelo, 'nls')) {
    
    datos_tmp <- eval(modelo$data)
    npar_mod <- length(coef(modelo))
    nro_datos <- nrow(datos_tmp)
    nro_obs_mod_i <- k.ini:nro_datos
    
    # La matriz para almacenar los valores de los coeficientes tiene en la
    # primera columna el número de observaciones usadas. En las siguientes
    # se almacenan los valores de los coeficientes con sus interv de conf
    mtrx_coef <- matrix(data = NA,
                        nrow = length(nro_obs_mod_i),
                        ncol = npar_mod * 3 + 1)
    
    idx_coef <- seq(2, ncol(mtrx_coef), 3)
    idx_lic <- seq(3, ncol(mtrx_coef), 3)
    idx_lsc <- seq(4, ncol(mtrx_coef), 3)
    
    colnames(mtrx_coef) <- paste('col', 1:ncol(mtrx_coef), sep = '.')
    colnames(mtrx_coef)[1] <- 'k'
    colnames(mtrx_coef)[idx_coef] <- names(coef(modelo))
    colnames(mtrx_coef)[idx_lic] <- paste(names(coef(modelo)), 'lic', sep = '.')
    colnames(mtrx_coef)[idx_lsc] <- paste(names(coef(modelo)), 'lsc', sep = '.')
    
    # lista de valores iniciales para modelo nls
    coef_aux_ini <- as.list(coef(modelo))
    
    for (i in nro_obs_mod_i) {
      mod_tmp <- nls(formula(modelo), start = coef_aux_ini,
                     data = datos_tmp[1:i, ])
      sum_tmp <- summary(mod_tmp)
      coef_tmp <- coef(mod_tmp)
      t_0 <- qt(0.025, sum_tmp$df[2])
      lic <- coef_tmp - t_0 * sum_tmp$coefficients[, 2]
      lsc <- coef_tmp + t_0 * sum_tmp$coefficients[, 2]
      mtrx_coef[i - k_ini + 1, 1] <- i
      mtrx_coef[i - k_ini + 1, idx_coef] <- coef_tmp
      mtrx_coef[i - k_ini + 1, idx_lic] <- lic
      mtrx_coef[i - k_ini + 1, idx_lsc] <- lsc
      
      coef_aux_ini <- as.list(coef(mod_tmp)) # valores iniciales próxima iter
    }
    
  } else {
    return('ERROR: Objeto no válido.')
  }
  
  return(mtrx_coef)
}


## Lectura de datos ------------------------------------------------------------

datos_orig <- fread('../datos/datos.csv', header = FALSE, skip = 11,
                         sep = ';', dec = ',',
                         colClasses = c(rep("NULL", 12),
                                        "numeric",
                                        rep("NULL", 5)))
setnames(datos_orig, 'ivand')
datos_orig[, log_ivand := log(ivand)]

# Definir algunos parámetros
nro_datos <- datos_orig[, .N]
n_lag <- 24 # número de rezagos en la ACF y PACF

datos_orig_ts <- ts(datos_orig$ivand, freq = 4, start = c(2000, 1))
datos_log_ts <- ts(datos_orig$log_ivand, freq = 4, start = c(2000, 1))

# organizamos los datos en un data frame y creamos otras variables de interés
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

# Gráficos de la serie observada
serie_orig <- ggplot(data = datos_orig, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand)) +
  labs(x = '',
       y = 'IVA no deducible [1000\'s mill pesos]') +
  scale_y_continuous(labels = dollar)
serie_orig

# Gráfico del logaritmo de la serie
serie_log <- ggplot(data = datos_orig, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = log_ivand)) +
  labs(x = '',
       y = 'log(IVA no deducible [1000\'s mill pesos])')
serie_log

# Guardar los dos gráficos
pdf(file = '../img/series.pdf', width = 7, height = 3.5)

grid.arrange(serie_orig, serie_log, nrow = 1, ncol = 2)

dev.off()

# Guardar gráfico de la ACF
pdf(file = '../img/acf_log_serie.pdf', width = 6, height = 4)

opar <- par(no.readonly = TRUE)
par(mar = c(4, 4, 1, 1))

# ACF del logaritmo de la serie
datos_orig[, stats::acf(log_ivand, lag.max = n_lag, ci.type = 'ma',
                        ylim = c(-1, 1), xlim = c(1, n_lag), main = '',
                        xlab = 'Rezago (k)')]
par(opar)

dev.off()


## Ajuste de modelos globales --------------------------------------------------

# Datos para el ajuste con validación cruzada. Omitimos los 4 últimos trimestres
datos_ajuste <- datos_orig[1:(nro_datos - 4)]

# Definir variables indicadoras para la estacionalidad
datos_ajuste[, c('I2', 'I3', 'I4') := list(
  ifelse(tri == 'T2', 1, 0),
  ifelse(tri == 'T3', 1, 0),
  ifelse(tri == 'T4', 1, 0)
)]

# Ajuste modelo log-cuadrático estacional
# $log(Y_t)=\beta_0+\beta_1t+\beta_2t^2+\delta_2I_{2t}+\delta_3I_{3t}+\delta_4I_{4t}+E_t$
modelo <- lm(log_ivand ~ t + I(t ^ 2) + I2 + I3 + I4, data = datos_ajuste)
summary(modelo)

# Imprimir tabla en formato LaTeX
formato_tab_latex <- c('s', 'fg', 'fg', 'g', 'g')
xtable(summary(modelo),
       caption = 'Ajuste modelo log-cuadrático estacional validación cruzada.',
       label = 'tab:modelo_aju',
       align = 'lrrrr',
       display = formato_tab_latex)

# $\sqrt{MSE}$ estimación de $\sigma$
sigma_mod <- summary(modelo)$sigma

# factor de corrección lognormal $exp(0.5\sigma^2)$
(factor_corr <- exp(0.5 * sigma_mod ^ 2))

datos_ajuste[, aju := fitted(modelo)] # valores ajustados escala logarítmica
datos_ajuste[, res := residuals(modelo)] # error $E_t$ escala logarítmica

datos_ajuste[, aju_orig := exp(aju) * factor_corr] # ajustados escala original
datos_ajuste[, res_orig := ivand - aju_orig] # residuales escala original

# Criterios de información $exp(C_n^*(p))$
exp(C_p(datos_ajuste[, res_orig], type = 'AIC', p = 6))
exp(C_p(datos_ajuste[, res_orig], type = 'BIC', p = 6))

exp(coef(modelo) * factor_corr) # parámetros escala original

# valores críticos $t_0$
qt(0.01 / 2, 49) # $\alpha = 1 \%$
qt(0.05 / 2, 49) # $\alpha = 5 \%$
qt(0.10 / 2, 49) # $\alpha = 10 \%$
qt(0.15 / 2, 49) # $\alpha = 15 \%$
qt(0.20 / 2, 49) # $\alpha = 20 \%$

# Gráfico de la serie  observada y la ajustada en escala original
orig_vs_ajuste <- ggplot(data = datos_ajuste, aes(x = fecha)) + theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Guardar gráfico
pdf(file = '../img/orig_vs_ajuste.ped', width = 5, height = 5)

orig_vs_ajuste

dev.off()

# Pronósticos para validación cruzada
nro_datos_aju <- datos_ajuste[, .N]
t_pron <- (nro_datos_aju + 1):nro_datos
I4_pron <- c(1, 0, 0, 0)
I2_pron <- c(0, 0, 1, 0)
I3_pron <- c(0, 0, 0, 1)
pred_mod <- predict(modelo,
                    newdata = data.frame(t = t_pron, I2 = I2_pron,
                                         I3 = I3_pron, I4 = I4_pron),
                    interval = 'prediction', level = 0.95)

pred_mod_esc_orig <- exp(pred_mod) * factor_corr

# Amplitud promedio IC
mean(pred_mod_esc_orig[, 3] - pred_mod_esc_orig[, 2])

# Medidas de precisión de los pronósticos
accuracy(datos_orig$pron_p[t_pron], datos_orig$ivand[t_pron])

datos_orig[, pron_p := c(rep(NA, nro_datos_aju), pred_mod_esc_orig[, 1])]
datos_orig[, pron_ici := c(rep(NA, nro_datos_aju), pred_mod_esc_orig[, 2])]
datos_orig[, pron_ics := c(rep(NA, nro_datos_aju), pred_mod_esc_orig[, 3])]

# Gráfico de la serie observada, ajustada y pronósticos con IC 95 %
serie_y_prono <- ggplot(data = datos_orig, aes(x = fecha)) + theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = pron_p, colour = 'Pronóstico'), linetype = 2) +
  geom_line(aes(y = pron_ici), linetype = 3) +
  geom_line(aes(y = pron_ics), linetype = 3) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Pronóstico' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Guardar gráficos de series observadas, ajuste y pronóstico
pdf(file = '../img/series_aju_y_prono.pdf', width = 7, height = 3.5)

grid.arrange(orig_vs_ajuste, serie_y_prono, nrow = 1, ncol  = 2)

dev.off()


## Postulación de modeloa ARMA para el error E_t -------------------------------

sd_resid <- datos_ajuste[, sd(res)]

# Gráficos de los residuales
res_vs_t <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res)) +
  geom_hline(yintercept = c(-2 * sd_resid, 0, 2 * sd_resid),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju, y = res)) +
  geom_hline(yintercept = c(-2 * sd_resid, 0, 2 * sd_resid),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales.pdf', width = 7, height = 3.5)

grid.arrange(res_vs_t, res_vs_aju, nrow = 1, ncol = 2)

dev.off()

# ACF y PACF de los errores
pdf(file = '../img/acf_pacf_res.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res, lag.max = n_lag,
                   ylim = c(-1, 1), xlim = c(1, n_lag),
                   main = '', xlab = 'Rezago (k)', ylab = 'PACF')]

par(opar)

dev.off()

# Prueba Durbin-Watson
xtable(
  pruebaDW(modelo),
  caption = 'Resultados prueba de orden 1 Durbin-Watson.',
  label = 'tab:res_mod_db'
)

# Prueba Box-Pierce
xtable(
  pruebaBPJB(datos_ajuste[, res], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box
xtable(
  pruebaBPJB(datos_ajuste[, res], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Método EACF
res_ts <- ts(datos_ajuste[, res], freq = 4, start = c(2000, 1))
eacf(res_ts, ar.max = 12, ma.max = 12)

# Método armasubsets()
pdf(file = '../img/armasubset.pdf', width = 7, height = 4)
plot(armasubsets(datos_ajuste[, res], y.name = 'Et',
                 nar = 12, nma = 12, ar.method = 'ml'))
dev.off()

# Método SelectModel()
SelectModel(datos_ajuste[, res],
            lag.max = n_lag, Criterion = 'AIC', ARModel = 'AR')
SelectModel(datos_ajuste[, res],
            lag.max = n_lag, Criterion = 'BIC', ARModel = 'AR')

# Método auto.arima()
auto.arima(datos_ajuste[, res], max.p = n_lag, max.q = n_lag, ic = 'aic')
auto.arima(res_ts, max.p = n_lag, max.q = n_lag, ic = 'aic')
auto.arima(datos_ajuste[, res], max.p = n_lag, max.q = n_lag, ic = 'bic')
auto.arima(res_ts, max.p = n_lag, max.q = n_lag, ic = 'bic')


## Ajuste de modelos ARIMA con validación cruzada ------------------------------

# Definir matriz de variables regresoras
mat_reg <- matrix(cbind(
  t = datos_ajuste[, t],
  t2 = datos_ajuste[, t ^ 2],
  I2 = datos_ajuste[, I2],
  I3 = datos_ajuste[, I3],
  I4 = datos_ajuste[, I4]
), ncol = 5)


# AR(1) ----
mod_ar1 <- Arima(datos_log_ts[1:nro_datos_aju],
                 order = c(1, 0, 0), xreg = mat_reg, method = 'ML')
summary(mod_ar1)

resumen_arima(mod_ar1)

xtable(resumen_arima(mod_ar1),
       caption = 'Resumen modelo con errores estructurales AR(1).',
       label = 'tab:mod_ar1',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar1 := resid(mod_ar1)]
datos_ajuste[, aju_ar1 := fitted(mod_ar1)]

factor_corr_ar1 <- exp(0.5 * mod_ar1$sigma2)

datos_ajuste[, aju_orig_ar1 := exp(aju_ar1) * factor_corr_ar1]
datos_ajuste[, res_orig_ar1 := ivand - aju_orig_ar1]

sd_resid_ar1 <- datos_ajuste[, sd(res_ar1)]

exp(C_p(datos_ajuste[, res_orig_ar1], type = 'AIC', p = 7))
exp(C_p(datos_ajuste[, res_orig_ar1], type = 'BIC', p = 7))

# Gráfico de la serie observada y ajustada
orig_vs_ajuste_ar1 <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar1, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo  AR(1)
res_vs_t_ar1 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar1)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar1, 0, 2 * sd_resid_ar1),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar1 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar1, y = res_ar1)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar1, 0, 2 * sd_resid_ar1),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar1.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar1, ggqqplot(mod_ar1),
             res_vs_t_ar1, res_vs_aju_ar1,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo AR(1)
pdf(file = '../img/acf_pacf_res_ar1.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar1, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar1, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce AR(1)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar1], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box AR(1)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar1], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos AR(1)
t_L <- (nro_datos_aju + 1):nro_datos
xreg_pron <- cbind(t = t_L, t2 = t_L ^ 2,
                   I2 = I2_pron, I3 = I3_pron, I4 = I4_pron)

pred_mod_ar1 <- forecast(mod_ar1, xreg = xreg_pron, level = 0.95)
pred_mod_ar1 <- exp(as.data.frame(pred_mod_ar1)) * factor_corr_ar1
pred_mod_ar1 <- ts(pred_mod_ar1, freq = 4, start = c(2000, 1))

# Precisión pronósticos AR(1)
accuracy(pred_mod_ar1, datos_orig_ts[t_L])

# Amplitud promedio pronósticos AR(1)
mean(pred_mod_ar1[, 3] - pred_mod_ar1[, 2])


# AR(5) ----
mod_ar5 <- Arima(datos_log_ts[1:nro_datos_aju],
             order = c(5, 0, 0), xreg = mat_reg, method = 'ML',
             fixed = c(NA, 0, 0, 0, NA, # ar
                       NA, NA, NA, NA, NA, NA)) # reg
summary(mod_ar5)

resumen_arima(mod_ar5)

xtable(resumen_arima(mod_ar5),
       caption = 'Resumen modelo con errores estructurales AR(5).',
       label = 'tab:mod_ar5',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar5 := resid(mod_ar5)]
datos_ajuste[, aju_ar5 := fitted(mod_ar5)]

factor_corr_ar5 <- exp(0.5 * mod_ar5$sigma2)

datos_ajuste[, aju_orig_ar5 := exp(aju_ar5) * factor_corr_ar5]
datos_ajuste[, res_orig_ar5 := ivand - aju_orig_ar5]

(coef_orig_ar5 <- exp(mod_ar5$coef) * factor_corr_ar5)

sd_resid_ar5 <- datos_ajuste[, sd(res_ar5)]

exp(C_p(datos_ajuste[, res_orig_ar5], type = 'AIC', p = 8))
exp(C_p(datos_ajuste[, res_orig_ar5], type = 'BIC', p = 8))

# Gráfico de la serie observada y ajustada AR(5)
orig_vs_ajuste_ar5 <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar5, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo  AR(5)
res_vs_t_ar5 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar5)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar5, 0, 2 * sd_resid_ar5),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar5 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar5, y = res_ar5)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar5, 0, 2 * sd_resid_ar5),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar5.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar5, ggqqplot(mod_ar5),
             res_vs_t_ar5, res_vs_aju_ar5,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo AR(5)
pdf(file = '../img/acf_pacf_res_ar5.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar5, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar5, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce AR(5)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar5], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box AR(5)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar5], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos AR(5)
pred_mod_ar5 <- forecast(mod_ar5, xreg = xreg_pron, level = 0.95)
pred_mod_ar5 <- exp(as.data.frame(pred_mod_ar5)) * factor_corr_ar5
pred_mod_ar5 <- ts(pred_mod_ar5, freq = 4, start = c(2000, 1))

# Precisión pronósticos AR(5)
accuracy(pred_mod_ar5, datos_orig_ts[t_L])

# Amplitud promedio pronósticos AR(5)
mean(pred_mod_ar5[, 3] - pred_mod_ar5[, 2])


# ARMA(5,6) ----
mod_ar5ma6 <- Arima(datos_log_ts[1:nro_datos_aju],
                    order = c(5, 0, 6), xreg = mat_reg, method = 'ML',
                    fixed = c(NA, 0, 0, 0, NA, # ar
                              0, 0, 0, 0, 0, NA, # ma
                              NA, NA, NA, NA, NA, NA)) # reg
summary(mod_ar5ma6)

resumen_arima(mod_ar5ma6)[c(1, 5, 11:17), ]

xtable(resumen_arima(mod_ar5ma6)[c(1, 5, 11:17), ],
       caption = 'Resumen modelo con errores estructurales ARMA(5,6).',
       label = 'tab:mod_ar5ma6',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar5ma6 := resid(mod_ar5ma6)]
datos_ajuste[, aju_ar5ma6 := fitted(mod_ar5ma6)]

factor_corr_ar5ma6 <- exp(0.5 * mod_ar5ma6$sigma2)

datos_ajuste[, aju_orig_ar5ma6 := exp(aju_ar5ma6) * factor_corr_ar5ma6]
datos_ajuste[, res_orig_ar5ma6 := ivand - aju_orig_ar5ma6]

sd_resid_ar5ma6 <- datos_ajuste[, sd(res_ar5ma6)]

exp(C_p(datos_ajuste[, res_orig_ar5ma6], type = 'AIC', p = 9))
exp(C_p(datos_ajuste[, res_orig_ar5ma6], type = 'BIC', p = 9))

# Gráfico de la serie observada y ajustada ARMA(5,6)
orig_vs_ajuste_ar5ma6 <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar5ma6, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo  ARMA(5,6)
res_vs_t_ar5ma6 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar5ma6)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar5ma6, 0, 2 * sd_resid_ar5ma6),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar5ma6 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar5ma6, y = res_ar5ma6)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar5ma6, 0, 2 * sd_resid_ar5ma6),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar5ma6.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar5ma6, ggqqplot(mod_ar5ma6),
             res_vs_t_ar5ma6, res_vs_aju_ar5ma6,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo ARMA(5,6)
pdf(file = '../img/acf_pacf_res_ar5ma6.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar5ma6, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar5ma6, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce ARMA(5,6)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar5ma6], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box ARMA(5,6)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar5ma6], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos ARMA(5,6)
pred_mod_ar5ma6 <- forecast(mod_ar5ma6, xreg = xreg_pron, level = 0.95)
pred_mod_ar5ma6 <- exp(as.data.frame(pred_mod_ar5ma6)) * factor_corr_ar5ma6
pred_mod_ar5ma6 <- ts(pred_mod_ar5ma6, freq = 4, start = c(2000, 1))

# Precisión pronósticos ARMA(5,6)
accuracy(pred_mod_ar5ma6, datos_orig_ts[t_L])

# Amplitud promedio pronósticos ARMA(5,6)
mean(pred_mod_ar5ma6[, 3] - pred_mod_ar5ma6[, 2])


# AR(3)MA(1)[4] ----

mod_ar3ma1s <- Arima(datos_log_ts[1:nro_datos_aju],
                     order = c(3, 0, 0),
                     seasonal = list(order = c(0, 0, 1), period = 4),
                     xreg = mat_reg, method = 'ML',
                     include.mean = FALSE)
summary(mod_ar3ma1s)

resumen_arima(mod_ar3ma1s)

xtable(resumen_arima(mod_ar3ma1s)[c(1, 5, 11:17), ],
       caption = 'Resumen modelo con errores estructurales ARMA(5,6).',
       label = 'tab:mod_ar3ma1s',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar3ma1s := resid(mod_ar3ma1s)]
datos_ajuste[, aju_ar3ma1s := fitted(mod_ar3ma1s)]

factor_corr_ar3ma1s <- exp(0.5 * mod_ar3ma1s$sigma2)

datos_ajuste[, aju_orig_ar3ma1s := exp(aju_ar3ma1s) * factor_corr_ar3ma1s]
datos_ajuste[, res_orig_ar3ma1s := ivand - aju_orig_ar3ma1s]

sd_resid_ar3ma1s <- datos_ajuste[, sd(res_ar3ma1s)]

exp(C_p(datos_ajuste[, res_orig_ar3ma1s], type = 'AIC', p = 9))
exp(C_p(datos_ajuste[, res_orig_ar3ma1s], type = 'BIC', p = 9))

# Gráfico de la serie observada y ajustada AR(3)MA(1)[4]
orig_vs_ajuste_ar3ma1s <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar3ma1s, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo AR(3)MA(1)[4]
res_vs_t_ar3ma1s <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar3ma1s)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar3ma1s, 0, 2 * sd_resid_ar3ma1s),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar3ma1s <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar3ma1s, y = res_ar3ma1s)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar3ma1s, 0, 2 * sd_resid_ar3ma1s),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar3ma1s.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar3ma1s, ggqqplot(mod_ar3ma1s),
             res_vs_t_ar3ma1s, res_vs_aju_ar3ma1s,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo AR(3)MA(1)[4]
pdf(file = '../img/acf_pacf_res_ar3ma1s.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar3ma1s, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar3ma1s, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce AR(3)MA(1)[4]
xtable(
  pruebaBPJB(datos_ajuste[, res_ar3ma1s], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box AR(3)MA(1)[4]
xtable(
  pruebaBPJB(datos_ajuste[, res_ar3ma1s], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos AR(3)MA(1)[4]
pred_mod_ar3ma1s <- forecast(mod_ar3ma1s, xreg = xreg_pron, level = 0.95)
pred_mod_ar3ma1s <- exp(as.data.frame(pred_mod_ar3ma1s)) * factor_corr_ar3ma1s
pred_mod_ar3ma1s <- ts(pred_mod_ar3ma1s, freq = 4, start = c(2000, 1))

# Precisión pronósticos AR(3)MA(1)[4]
accuracy(pred_mod_ar3ma1s, datos_orig_ts[t_L])

# Amplitud promedio pronósticos AR(3)MA(1)[4]
mean(pred_mod_ar3ma1s[, 3] - pred_mod_ar3ma1s[, 2])


# AR(2)MA(1)[4] ----

mod_ar2ma1s <- Arima(datos_log_ts[1:nro_datos_aju],
                     order = c(2, 0, 0),
                     seasonal = list(order = c(0, 0, 1), period = 4),
                     xreg = mat_reg, method = 'ML')

summary(mod_ar2ma1s)

resumen_arima(mod_ar2ma1s)

# res_ar2ma1s_x13 <- scan("")
# 0.004 -0.010 -0.015 0.001 -0.033 0.012 -0.015 0.003 -0.035 0.019 0.003 -0.011 0.026 -0.035 0.006 0.018 -0.040 0.006 0.002 0.011 -0.023 0.033 -0.011 -0.010 0.028 0.008 0.043 0.009 0.014 -0.011 0.006 0.014 -0.007 0.006 0.017 -0.025 -0.035 -0.034 -0.018 -0.014 0.018 0.003 0.016 0.009 0.014 0.018 -0.019 -0.001 0.015 -0.034 -0.010 -0.005 0.040 -0.002 -0.011
# 
# t_46 <- pt(c(159.88, 11.75, 3.86, 3.72, 4.32, 17.21,
#              1.1556/0.1204, 0.3857/0.12368, 0.6733/0.10680), 46,
#            lower.tail = FALSE)
# 
# datos_ajuste[, res_ar2ma1s_x13 := res_ar2ma1s_x13]

xtable(resumen_arima(mod_ar2ma1s),
       caption = 'Resumen modelo con errores estructurales ARMA(5,6).',
       label = 'tab:mod_ar2ma1s',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar2ma1s := resid(mod_ar2ma1s)]
datos_ajuste[, aju_ar2ma1s := fitted(mod_ar2ma1s)]

factor_corr_ar2ma1s <- exp(0.5 * mod_ar2ma1s$sigma2)

datos_ajuste[, aju_orig_ar2ma1s := exp(aju_ar2ma1s) * factor_corr_ar2ma1s]
datos_ajuste[, res_orig_ar2ma1s := ivand - aju_orig_ar2ma1s]

exp(mod_ar2ma1s$coef) * factor_corr_ar2ma1s

sd_resid_ar2ma1s <- datos_ajuste[, sd(res_ar2ma1s)]

exp(C_p(datos_ajuste[, res_orig_ar2ma1s], type = 'AIC', p = 9))
exp(C_p(datos_ajuste[, res_orig_ar2ma1s], type = 'BIC', p = 9))

# Gráfico de la serie observada y ajustada AR(2)MA(1)[4]
orig_vs_ajuste_ar2ma1s <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar2ma1s, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo AR(2)MA(1)[4]
res_vs_t_ar2ma1s <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar2ma1s)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar2ma1s, 0, 2 * sd_resid_ar2ma1s),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar2ma1s <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar2ma1s, y = res_ar2ma1s)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar2ma1s, 0, 2 * sd_resid_ar2ma1s),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar2ma1s.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar2ma1s, ggqqplot(mod_ar2ma1s),
             res_vs_t_ar2ma1s, res_vs_aju_ar2ma1s,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo AR(2)MA(1)[4]
pdf(file = '../img/acf_pacf_res_ar2ma1s.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar2ma1s, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar2ma1s, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce AR(2)MA(1)[4]
xtable(
  pruebaBPJB(datos_ajuste[, res_ar2ma1s], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box AR(2)MA(1)[4]
xtable(
  pruebaBPJB(datos_ajuste[, res_ar2ma1s], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos AR(2)MA(1)[4]
pred_mod_ar2ma1s <- forecast(mod_ar2ma1s, xreg = xreg_pron, level = 0.95)
pred_mod_ar2ma1s <- exp(as.data.frame(pred_mod_ar2ma1s)) * factor_corr_ar2ma1s
pred_mod_ar2ma1s <- ts(pred_mod_ar2ma1s, freq = 4, start = c(2000, 1))

# Precisión pronósticos AR(2)MA(1)[4]
accuracy(pred_mod_ar2ma1s, datos_orig_ts[t_L])

# Amplitud promedio pronósticos AR(2)MA(1)[4]
mean(pred_mod_ar2ma1s[, 3] - pred_mod_ar2ma1s[, 2])


# AR(4) ----

mod_ar4 <- Arima(datos_log_ts[1:nro_datos_aju],
                     order = c(4, 0, 0),
                     xreg = mat_reg, method = 'ML')
summary(mod_ar4)

resumen_arima(mod_ar4)

xtable(resumen_arima(mod_ar4)[c(1, 5, 11:17), ],
       caption = 'Resumen modelo con errores estructurales ARMA(5,6).',
       label = 'tab:mod_ar4',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar4 := resid(mod_ar4)]
datos_ajuste[, aju_ar4 := fitted(mod_ar4)]

factor_corr_ar4 <- exp(0.5 * mod_ar4$sigma2)

datos_ajuste[, aju_orig_ar4 := exp(aju_ar4) * factor_corr_ar4]
datos_ajuste[, res_orig_ar4 := ivand - aju_orig_ar4]

sd_resid_ar4 <- datos_ajuste[, sd(res_ar4)]

exp(C_p(datos_ajuste[, res_orig_ar4], type = 'AIC', p = 9))
exp(C_p(datos_ajuste[, res_orig_ar4], type = 'BIC', p = 9))

# Gráfico de la serie observada y ajustada AR(4)
orig_vs_ajuste_ar4 <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar4, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo AR(4)
res_vs_t_ar4 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar4)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar4, 0, 2 * sd_resid_ar4),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar4 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar4, y = res_ar4)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar4, 0, 2 * sd_resid_ar4),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar4.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar4, ggqqplot(mod_ar4),
             res_vs_t_ar4, res_vs_aju_ar4,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo AR(4)
pdf(file = '../img/acf_pacf_res_ar4.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar4, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar4, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce AR(4)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar4], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box AR(4)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar4], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos AR(4)
pred_mod_ar4 <- forecast(mod_ar4, xreg = xreg_pron, level = 0.95)
pred_mod_ar4 <- exp(as.data.frame(pred_mod_ar4)) * factor_corr_ar4
pred_mod_ar4 <- ts(pred_mod_ar4, freq = 4, start = c(2000, 1))

# Precisión pronósticos AR(4)
accuracy(pred_mod_ar4, datos_orig_ts[t_L])

# Amplitud promedio pronósticos AR(4)
mean(pred_mod_ar4[, 3] - pred_mod_ar4[, 2])


# AR(5) full ----

mod_ar5f <- Arima(datos_log_ts[1:nro_datos_aju],
                 order = c(5, 0, 0),
                 xreg = mat_reg, method = 'ML')
summary(mod_ar5f)

resumen_arima(mod_ar5f)

pt(c(0.9352/0.1240, 0.0252/0.1821, 6.49, 5.19), 44, lower.tail = FALSE)

xtable(resumen_arima(mod_ar5f)[c(1, 5, 11:17), ],
       caption = 'Resumen modelo con errores estructurales ARMA(5,6).',
       label = 'tab:mod_ar5f',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar5f := resid(mod_ar5f)]
datos_ajuste[, aju_ar5f := fitted(mod_ar5f)]

factor_corr_ar5f <- exp(0.5 * mod_ar5f$sigma2)

datos_ajuste[, aju_orig_ar5f := exp(aju_ar5f) * factor_corr_ar5f]
datos_ajuste[, res_orig_ar5f := ivand - aju_orig_ar5f]

sd_resid_ar5f <- datos_ajuste[, sd(res_ar5f)]

exp(C_p(datos_ajuste[, res_orig_ar5f], type = 'AIC', p = 9))
exp(C_p(datos_ajuste[, res_orig_ar5f], type = 'BIC', p = 9))

# Gráfico de la serie observada y ajustada AR(5) full
orig_vs_ajuste_ar5f <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar5f, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo AR(5) full
res_vs_t_ar5f <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar5f)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar5f, 0, 2 * sd_resid_ar5f),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar5f <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar5f, y = res_ar5f)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar5f, 0, 2 * sd_resid_ar5f),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar5f.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar5f, ggqqplot(mod_ar5f),
             res_vs_t_ar5f, res_vs_aju_ar5f,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo AR(5) full
pdf(file = '../img/acf_pacf_res_ar5f.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar5f, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar5f, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce AR(5) full
xtable(
  pruebaBPJB(datos_ajuste[, res_ar5f], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box AR(5) full
xtable(
  pruebaBPJB(datos_ajuste[, res_ar5f], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos AR(5) full
pred_mod_ar5f <- forecast(mod_ar5f, xreg = xreg_pron, level = 0.95)
pred_mod_ar5f <- exp(as.data.frame(pred_mod_ar5f)) * factor_corr_ar5f
pred_mod_ar5f <- ts(pred_mod_ar5f, freq = 4, start = c(2013, 4))

# Precisión pronósticos AR(5) full
accuracy(pred_mod_ar5f, datos_orig_ts[t_L])

# Amplitud promedio pronósticos AR(5) full
mean(pred_mod_ar5f[, 3] - pred_mod_ar5f[, 2])


# ARMA(1,3) ----

mod_ar1ma3 <- Arima(datos_log_ts[1:nro_datos_aju],
                  order = c(1, 0, 3),
                  xreg = mat_reg, method = 'ML')
summary(mod_ar1ma3)

resumen_arima(mod_ar1ma3)

xtable(resumen_arima(mod_ar1ma3)[c(1, 5, 11:17), ],
       caption = 'Resumen modelo con errores estructurales ARMA(5,6).',
       label = 'tab:mod_ar1ma3',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar1ma3 := resid(mod_ar1ma3)]
datos_ajuste[, aju_ar1ma3 := fitted(mod_ar1ma3)]

factor_corr_ar1ma3 <- exp(0.5 * mod_ar1ma3$sigma2)

datos_ajuste[, aju_orig_ar1ma3 := exp(aju_ar1ma3) * factor_corr_ar1ma3]
datos_ajuste[, res_orig_ar1ma3 := ivand - aju_orig_ar1ma3]

sd_resid_ar1ma3 <- datos_ajuste[, sd(res_ar1ma3)]

exp(C_p(datos_ajuste[, res_orig_ar1ma3], type = 'AIC', p = 9))
exp(C_p(datos_ajuste[, res_orig_ar1ma3], type = 'BIC', p = 9))

# Gráfico de la serie observada y ajustada ARMA(1,3)
orig_vs_ajuste_ar1ma3 <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar1ma3, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo ARMA(1,3)
res_vs_t_ar1ma3 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar1ma3)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar1ma3, 0, 2 * sd_resid_ar1ma3),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar1ma3 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar1ma3, y = res_ar1ma3)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar1ma3, 0, 2 * sd_resid_ar1ma3),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar1ma3.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar1ma3, ggqqplot(mod_ar1ma3),
             res_vs_t_ar1ma3, res_vs_aju_ar1ma3,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo ARMA(1,3)
pdf(file = '../img/acf_pacf_res_ar1ma3.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar1ma3, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar1ma3, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce ARMA(1,3)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar1ma3], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box ARMA(1,3)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar1ma3], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos ARMA(1,3)
pred_mod_ar1ma3 <- forecast(mod_ar1ma3, xreg = xreg_pron, level = 0.95)
pred_mod_ar1ma3 <- exp(as.data.frame(pred_mod_ar1ma3)) * factor_corr_ar1ma3
pred_mod_ar1ma3 <- ts(pred_mod_ar1ma3, freq = 4, start = c(2000, 1))

# Precisión pronósticos ARMA(1,3)
accuracy(pred_mod_ar1ma3, datos_orig_ts[t_L])

# Amplitud promedio pronósticos ARMA(1,3)
mean(pred_mod_ar1ma3[, 3] - pred_mod_ar1ma3[, 2])


# AR(2) ----

mod_ar2 <- Arima(datos_log_ts[1:nro_datos_aju],
                    order = c(2, 0, 0),
                    xreg = mat_reg, method = 'ML')
summary(mod_ar2)

resumen_arima(mod_ar2)

xtable(resumen_arima(mod_ar2)[c(1, 5, 11:17), ],
       caption = 'Resumen modelo con errores estructurales ARMA(5,6).',
       label = 'tab:mod_ar2',
       display = c('s', 'fg', 'fg', 'g'))

datos_ajuste[, res_ar2 := resid(mod_ar2)]
datos_ajuste[, aju_ar2 := fitted(mod_ar2)]

factor_corr_ar2 <- exp(0.5 * mod_ar2$sigma2)

datos_ajuste[, aju_orig_ar2 := exp(aju_ar2) * factor_corr_ar2]
datos_ajuste[, res_orig_ar2 := ivand - aju_orig_ar2]

sd_resid_ar2 <- datos_ajuste[, sd(res_ar2)]

exp(C_p(datos_ajuste[, res_orig_ar2], type = 'AIC', p = 9))
exp(C_p(datos_ajuste[, res_orig_ar2], type = 'BIC', p = 9))

# Gráfico de la serie observada y ajustada AR(2)
orig_vs_ajuste_ar2 <- ggplot(data = datos_ajuste, aes(x = fecha)) +
  theme_bw(10) +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = aju_orig_ar2, colour = 'Ajuste')) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Ajuste' = 'red')) +
  theme(legend.position = c(0.20, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)

# Gráficos de los residuales modelo AR(2)
res_vs_t_ar2 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_line(aes(x = fecha, y = res_ar2)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar2, 0, 2 * sd_resid_ar2),
             colour = 'blue', linetype = 2) +
  labs(x = 'Tiempo', y = expression(hat(symbol(E))[t]))

res_vs_aju_ar2 <- ggplot(data = datos_ajuste) + theme_bw(10) +
  geom_point(aes(x = aju_ar2, y = res_ar2)) +
  geom_hline(yintercept = c(-2 * sd_resid_ar2, 0, 2 * sd_resid_ar2),
             colour = 'blue', linetype = 2) +
  labs(x = 'Ajustados', y = expression(hat(symbol(E))[t]))

pdf(file = '../img/residuales_ar2.pdf', width = 7, height = 6)

grid.arrange(orig_vs_ajuste_ar2, ggqqplot(mod_ar2),
             res_vs_t_ar2, res_vs_aju_ar2,
             nrow = 2, ncol = 2)

dev.off()

# ACF y PACF de los errores modelo AR(2)
pdf(file = '../img/acf_pacf_res_ar2.pdf', width = 9, height = 3.5)

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)

datos_ajuste[, stats::acf(res_ar2, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = 'Rezago (k)')]

datos_ajuste[, pacf(res_ar2, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = 'Rezago (k)', ylab = 'PACF')]
par(opar)

dev.off()

# Prueba Box-Pierce AR(2)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar2], maxlag = n_lag, type = 'Box'),
  caption = 'Resultados prueba Box-Pierce.',
  label = 'tab:res_mod_bp'
)
# Prueba Ljung-Box AR(2)
xtable(
  pruebaBPJB(datos_ajuste[, res_ar2], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Pronósticos AR(2)
pred_mod_ar2 <- forecast(mod_ar2, xreg = xreg_pron, level = 0.95)
pred_mod_ar2 <- exp(as.data.frame(pred_mod_ar2)) * factor_corr_ar2
pred_mod_ar2 <- ts(pred_mod_ar2, freq = 4, start = c(2000, 1))

# Precisión pronósticos AR(2)
accuracy(pred_mod_ar2, datos_orig_ts[t_L])

# Amplitud promedio pronósticos AR(2)
mean(pred_mod_ar2[, 3] - pred_mod_ar2[, 2])

# Gráficos de las series de los tres modelos ajustados

pdf('../img/series_vs_ajuste.pdf', width = 10, height = 3.5)

grid.arrange(orig_vs_ajuste_ar5 + labs(title = 'Modelo 1 AR(5) restringido'),
             orig_vs_ajuste_ar2ma1s + labs(title = 'Modelo 2 AR(2)MA(1)[4]'),
             orig_vs_ajuste_ar5f + labs(title = 'Modelo 3 AR(5) completo'),
             nrow = 1, ncol = 3)

dev.off()

## Análisis de residuales ------------------------------------------------------

xtable(
  pruebaBPJB(datos_ajuste[, res_ar5], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

xtable(
  pruebaBPJB(datos_ajuste[, res_ar2ma1s], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

xtable(
  pruebaBPJB(datos_ajuste[, res_ar5f], maxlag = n_lag, type = 'Ljung'),
  caption = 'Resultados prueba Ljung-Box.',
  label = 'tab:res_mod_lb'
)

# Gráficos análisis de residuales de los tres modelos seleccionados
pdf('../img/residuales_def.pdf', width = 10, height = 10)

grid.arrange(res_vs_t_ar5 + labs(title = 'Modelo 1'),
             res_vs_aju_ar5 + labs(title = 'Modelo 1'),
             ggqqplot(mod_ar5) + labs(title = 'Modelo 1'),
             res_vs_t_ar2ma1s + labs(title = 'Modelo 2'),
             res_vs_aju_ar2ma1s + labs(title = 'Modelo 2'),
             ggqqplot(mod_ar2ma1s) + labs(title = 'Modelo 2'),
             res_vs_t_ar5f + labs(title = 'Modelo 3'),
             res_vs_aju_ar5f + labs(title = 'Modelo 3'),
             ggqqplot(mod_ar5f) + labs(title = 'Modelo 3'),
             nrow = 3, ncol = 3)

dev.off()

# Gráficos ACF y PACF de los tres modelos seleccionados
pdf('../img/acf_pacf_residuales_def.pdf', width = 10, height = 10)

opar <- par(no.readonly = TRUE)
par(mfrow = c(3, 2), mar = c(4, 4, 1, 1),
    cex = 0.8, cex.lab = 0.8, cex.axis = 0.8)


datos_ajuste[, stats::acf(res_ar5, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          xlab = '', ylab = '', main = '')]
title(main = '', xlab = 'Rezago (k)', ylab = 'ACF', line = 2)
mtext('Modelo 1', cex = 0.7)

datos_ajuste[, pacf(res_ar5, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = '', ylab = '')]
title(main = '', xlab = 'Rezago (k)', ylab = 'PACF', line = 2)
mtext('Modelo 1', cex = 0.7)

datos_ajuste[, stats::acf(res_ar2ma1s, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = '', ylab = '')]
title(main = '', xlab = 'Rezago (k)', ylab = 'ACF', line = 2)
mtext('Modelo 2', cex = 0.7)

datos_ajuste[, pacf(res_ar2ma1s, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = '', ylab = '')]
title(main = '', xlab = 'Rezago (k)', ylab = 'PACF', line = 2)
mtext('Modelo 2', cex = 0.7)

datos_ajuste[, stats::acf(res_ar5f, lag.max = n_lag, ci.type = 'ma',
                          ylim = c(-1, 1), xlim = c(1, n_lag),
                          main = '', xlab = '', ylab = '')]
title(main = '', xlab = 'Rezago (k)', ylab = 'ACF', line = 2)
mtext('Modelo 3', cex = 0.7)

datos_ajuste[, pacf(res_ar5f, lag.max = n_lag,
                    ylim = c(-1, 1), xlim = c(1, n_lag),
                    main = '', xlab = '', ylab = '')]
title(main = '', xlab = 'Rezago (k)', ylab = 'PACF', line = 2)
mtext('Modelo 3', cex = 0.7)

par(opar)

dev.off()


## Predicciones

dat_obs <- datos_orig[(nro_datos_aju + 1):.N, ivand]
pred_ar5 <- as.vector(pred_mod_ar5[,1])
pred_ar2ma1s <- as.vector(pred_mod_ar2ma1s[,1])
pred_ar5f <- as.vector(pred_mod_ar5f[,3])

plot(dat_obs, type = 'l')
lines(pred_ar5, col = 'red', lty = 1, pch = 5)
lines(pred_ar2ma1s, col = 'blue', lty = 1, pch = 5)
lines(pred_ar5f, col = 'green', lty = 1, pch = 5)

datos_pred <- data.frame(Yt = c(dat_obs, pred_ar5, pred_ar2ma1s, pred_ar5f),
                         Serie = c('Real', 'Real', 'Real', 'Real',
                                   'Mod 1', 'Mod 1', 'Mod 1', 'Mod 1',
                                   'Mod 2', 'Mod 2', 'Mod 2', 'Mod 2',
                                   'Mod 3', 'Mod 3', 'Mod 3', 'Mod 3'))

pdf('../img/graf_pronost.pdf', width = 6, height = 5)

ggplot(data = datos_pred) + theme_bw(12) +
  geom_line(aes(x = rep(1:4, 4), y = Yt, linetype = Serie)) +
  geom_point(aes(x = rep(1:4, 4), y = Yt, shape = Serie), size = 3) +
  labs(x = '', y = 'Iva no deducible [1000 de millones]')

dev.off()
