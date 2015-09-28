#------------------------------------------------------------------------------#
# Código R trabajo 1: Análisis IVA No Deducible en Colombia                    #
# Autor: Juan Pablo Calle Quintero                                             #
# Fecha: 24/09/2015                                                            #
# Correo: jpcalleq@unal.edu.co                                                 #
#------------------------------------------------------------------------------#

## Librerías requeridas --------------------------------------------------------
library(ggplot2) # gráficos de calidad
library(grid) # para poder utilizar la función unit() en ggplot2
library(gridExtra) # funciones auxiliares graficación múltiple, grid.arrange()
library(dplyr) # manipulación de datos
library(xtable) # exportar tablas en formato LaTeX
library(scales) # formatear los valores de los ejes de las gráficas.
library(forecast) # funciones de pronóstico para series de tiempo y mod lineales
library(strucchange) # análisis de estabilidad
# theme_set(theme_bw(base_size = 7)) # definir tamaño de fuente gráficos ggplot2

## Funciones auxiliares --------------------------------------------------------
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

# Gráfico en blanco para usar en función grid.arrange()
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

## Lectura de datos ------------------------------------------------------------
datos_orig <- read.table('../datos/datos.csv', header = FALSE, skip = 11,
                         sep=';', dec = ',',
                         colClasses = c(rep("NULL", 12),
                                        "numeric",
                                        rep("NULL", 5)))
names(datos_orig) <- 'ivand'
nro_datos <- nrow(datos_orig)
datos_orig_ts <- ts(datos_orig, freq = 4, start = c(2000, 1))
datos_log_ts <- ts(log(datos_orig), freq = 4, start = c(2000, 1))

# organizamos los datos y le damos el formato adecuado
fecha <- seq.Date(from = as.Date('2000-03-01', '%Y-%m-%d'),
                  length.out = 59,
                  by = 'quarter')
mes <- reorder(months.Date(fecha), c(rep(c(3, 6, 9, 12), 14), 3, 6, 9))
tri <- gl(n = 4, k = 1, length = 59, labels = c('T1', 'T2', 'T3', 'T4'))
log_ivand <- log(datos_orig)
names(log_ivand) <- 'log_ivand'
datos <- data.frame(mes, tri, fecha, datos_orig, log_ivand)
head(datos)


## Análisis exploratorio inicial -----------------------------------------------
summary(datos)
plot(datos_orig_ts)
abline(v = 2000:2015, lty = 2)
plot(datos_log_ts)
abline(v = 2000:2015, lty = 2)

## 1. Análisis descriptivo de la serie -----------------------------------------

# Gráficos de la serie
marcas_anios <- seq.Date(as.Date('2000-01-01'),
                         as.Date('2015-01-01'), by='year')
p1_serie_orig <- ggplot(data = datos, aes(x = fecha)) +
  theme_bw() +
  geom_line(aes(y = ivand)) +
  geom_vline(xintercept = fecha, linetype = 2, colour = 'red') +
  labs(x = '',
       y = 'IVA no deducible [1000\'s mill pesos]') +
  scale_y_continuous(labels = dollar)
p1_serie_orig

p1_serie_log <- ggplot(data = datos, aes(x = fecha)) +
  theme_bw() +
  geom_line(aes(y = log_ivand)) +
  labs(x = '',
       y = 'log(IVA no deducible [1000\'s mill pesos])')
p1_serie_log

grid.arrange(p1_serie_orig, p1_serie_log,
             layout_matrix = matrix(c(1, 2), nrow = 1, byrow = TRUE))

# Gráficos descomposición de la serie
ts_desc_orig <- decompose(datos_orig_ts, type = 'multiplicative')

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

# Gráfico de la tendencia multiplicativa
datos$trend_mult <- as.vector(ts_desc_orig$trend)
datos_trend <- datos[!is.na(datos$trend_mult), ]

p1_trend_mult <- ggplot(data = datos_trend, aes(x = fecha)) +
  theme_bw() +
  geom_line(aes(y = trend_mult)) +
  labs(x = '',
       y = expression(paste(S[t])))
p1_trend_mult

ts_desc_log <- decompose(datos_log_ts, type = 'additive')

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

# Gráfico de la tendencia aditiva  del logaritmo
datos$trend_add_log <- as.vector(ts_desc_log$trend)
datos_trend$trend_add_log <- datos[!is.na(datos$trend_add_log), 'trend_add_log']

p1_trend_add_log <- ggplot(data = datos_trend, aes(x = fecha)) +
  theme_bw() +
  geom_line(aes(y = trend_add_log)) +
  labs(x = '',
       y = expression(paste('log(', S[t], ')')))
p1_trend_add_log

# Box plot
prom_trim_orig <- summarise(group_by(datos, mes), prom = mean(ivand))
prom_trim_log <- summarise(group_by(datos, mes), prom = mean(log_ivand))

trimestres <- c('T1', 'T2', 'T3', 'T4')

p1_boxplot_trim_orig <- ggplot(data = datos) +
  theme_bw() +
  geom_boxplot(aes(x = mes, y = ivand)) +
  geom_line(data = prom_trim_orig, aes(x = 1:4, y = prom),
            size = 0.2, colour = 'red', linetype = 6) +
  scale_x_discrete(labels = trimestres) +
  labs(x = '', y = 'IVA no desc [1000\'s mill pesos]') +
  theme(
    axis.line = element_blank(),
    #panel.border = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_y_continuous(labels = dollar)

p1_boxplot_trim_orig

p1_boxplot_trim_log <- ggplot(data = datos) +
  theme_bw() +
  geom_boxplot(aes(x = mes, y = log_ivand)) +
  geom_line(data = prom_trim_log, aes(x = 1:4, y = prom),
            size = 0.2, colour = 'red', linetype = 6) +
  scale_x_discrete(labels = trimestres) +
  labs(x = '', y = 'log(IVA no desc [1000\'s mill pesos])') +
  theme(
    axis.line = element_blank(),
    #panel.border = element_blank(),
    axis.ticks = element_blank()
  )

p1_boxplot_trim_log

# Periodogramas
diff_ivand <- diff(datos$ivand)
diff_log_ivand <- diff(datos$log_ivand)

nro_datos_dif <- length(diff_ivand)

period_ivand <- abs(fft(diff_ivand) / sqrt(nro_datos_dif))^2
period_log_ivand <- abs(fft(diff_log_ivand) / sqrt(nro_datos_dif))^2

period_ivand_scale <- period_ivand * (4/nro_datos_dif)
period_log_ivand_scale <- period_log_ivand * (4/nro_datos_dif)

freq <- (0:(nro_datos_dif - 1)) / nro_datos_dif
plot(freq, period_ivand_scale, type = 'l')
plot(freq, period_log_ivand_scale, type = 'l')

periodogramas <- data.frame(freq,
                            ivand = period_ivand_scale,
                            log_ivand = period_log_ivand_scale)

p1_periodograma_ivand <- ggplot(data = periodogramas) +
  theme_bw() +
  geom_line(aes(x = freq, y = ivand)) +
  geom_vline(xintercept = (1:3)/4, linetype = 2, colour = 'red') +
  labs(y = 'Periodograma', x = 'Frecuencia')

p1_periodograma_ivand

p1_periodograma_log_ivand <- ggplot(data = periodogramas) +
  theme_bw() +
  geom_line(aes(x = freq, y = log_ivand)) +
  geom_vline(xintercept = (1:3)/4, linetype = 2, colour = 'red') +
  labs(y = 'Periodograma', x = 'Frecuencia')

p1_periodograma_log_ivand

grid.arrange(p1_trend_mult, p1_trend_add_log,
             p1_boxplot_trim_log, p1_periodograma_log_ivand,
             nrow = 2, ncol = 2)

## 2. Postulación de modelos ---------------------------------------------------
# modelo 1: logarítmico-cuadrático con estacionalidad
# modelo 2: exponencial-cuadrático con estacionalidad

## 3. Ajuste de modelos con validación cruzada ---------------------------------
t_ajuste <- 1:(nro_datos - 4) # dejamos los últimos 4 trimestres por fuera
datos_ajuste <- datos[t_ajuste, ]
datos_ajuste$t <- t_ajuste

# modelo 1: logarítmico-cúbico
mod1 <- lm(log_ivand ~ t + I(t^2) + tri, data = datos_ajuste)
summary(mod1)
xtable(summary(mod1))

resid_mod1 <- resid(mod1)
sigma_mod1 <- summary(mod1)$sigma
datos_ajuste$resid_mod1 <- resid_mod1
aju_mod1 <- fitted(mod1)
datos_ajuste$aju_mod1 <- aju_mod1
resid_mod1_esc_orig <- datos_ajuste$ivand - exp(aju_mod1) * exp(.5*sigma_mod1^2)

exp(C_p(resid_mod1_esc_orig, type = 'AIC', p = 6))
exp(C_p(resid_mod1_esc_orig, type = 'BIC', p = 6))

qt(0.025, 49) # valor crítico de t

# modelo 2: exponencial-cúbico
val_ini <- coef(mod1)

I2 <- ifelse(datos_ajuste$tri == 'T2', 1, 0)
I3 <- ifelse(datos_ajuste$tri == 'T3', 1, 0)
I4 <- ifelse(datos_ajuste$tri == 'T4', 1, 0)

datos_ajuste$I2 <- I2
datos_ajuste$I3 <- I3
datos_ajuste$I4 <- I4

mod2 <- nls(ivand ~ exp(b0 + b1*t + b2*I(t^2) + d2*I2 + d3*I3 + d4*I4),
            start = list(b0 = val_ini[1], b1 = val_ini[2], b2 = val_ini[3],
                         d2 = val_ini[4], d3 = val_ini[5], d4 = val_ini[6]),
            data = datos_ajuste)
summary(mod2)

resid_mod2 <- resid(mod2)
datos_ajuste$resid_mod2 <- resid_mod2
datos_ajuste$aju_mod2 <- fitted(mod2)
exp(C_p(resid_mod2, type = 'AIC', p = 6))
exp(C_p(resid_mod2, type = 'BIC', p = 6))

# modelo 1 sin segundo y tercer trimestre

mod1_sin23trim <- lm(log_ivand ~ t + I(t^2) + I4, data = datos_ajuste)
summary(mod1_sin23trim)
xtable(summary(mod1_sin23trim))

resid_mod1_sin23trim <- resid(mod1_sin23trim)
sigma_mod1_sin23trim <- summary(mod1_sin23trim)$sigma
datos_ajuste$resid_mod1_sin23trim <- resid_mod1_sin23trim
aju_mod1_sin23trim <- fitted(mod1_sin23trim)
datos_ajuste$aju_mod1_sin23trim <- aju_mod1_sin23trim
datos_ajuste$aju_mod1_esc_orig_sin23trim <-
  exp(aju_mod1_sin23trim) * exp(.5*sigma_mod1_sin23trim^2)
resid_mod1_esc_orig_sin23trim <-
  datos_ajuste$ivand - datos_ajuste$aju_mod1_esc_orig_sin23trim

# coeficientes en escala original
exp(coef(mod1_sin23trim)) * exp(0.5*sigma_mod1_sin23trim^2)

exp(C_p(resid_mod1_esc_orig_sin23trim, type = 'AIC', p = 6))
exp(C_p(resid_mod1_esc_orig_sin23trim, type = 'BIC', p = 6))

qt(0.025, 51) # valor crítico de t

# modelo 2 sin segundo y tercer trimestre
mod2_sin23trim <- nls(ivand ~ exp(b0 + b1*t + b2*I(t^2) + d4*I4),
                      start = list(b0 = val_ini[1], b1 = val_ini[2],
                                   b2 = val_ini[3], d4 = val_ini[6]),
                      data = datos_ajuste)
summary(mod2_sin23trim)

resid_mod2_sin23trim <- resid(mod2_sin23trim)
datos_ajuste$resid_mod2_sin23trim <- resid_mod2_sin23trim
datos_ajuste$aju_mod2_sin23trim <- fitted(mod2_sin23trim)
exp(C_p(resid_mod2_sin23trim, type = 'AIC', p = 6))
exp(C_p(resid_mod2_sin23trim, type = 'BIC', p = 6))

# modelo 3: Holt-Winters Multiplicativo


# modelo 4: LOESS transformado Yt / St_descomp
St <- ts_desc_orig$seasonal[1:4]
datos_ajuste$St <- rep(St, nro_datos - 4)[1:(nro_datos - 4)]
Yt_desc <- datos_ajuste$ivand / datos_ajuste$St


## 4. Análisis de residuales y validación de supuestos -------------------------

# análisis residuales modelo 1 sin trimestres 2 y 3
p4_base <- ggplot(data = datos_ajuste, aes(x = fecha)) + theme_bw()

p4.1_mod1 <- p4_base + geom_line(aes(y = resid_mod1_sin23trim)) +
  labs(x = 'Tiempo',
       y  = 'Residuales')
p4.1_mod1

p4.2_mod1 <- p4_base + geom_point(aes(x = aju_mod1_sin23trim,
                                     y = resid_mod1_sin23trim)) +
  labs(x = 'Ajustados',
       y = 'Residuales')
p4.2_mod1

p4.3_mod1 <- p4_base + geom_line(aes(y = ivand, colour = 'Real')) + 
  geom_line(aes(y = aju_mod1_esc_orig_sin23trim, colour = 'Ajuste')) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') +
  scale_colour_manual("",
                      values = c('Real'='black', 'Ajuste'='red')) +
  theme(legend.position = c(0.15, 0.85), legend.background = element_blank())
p4.3_mod1

grid.arrange(p4.3_mod1, blankPlot, p4.1_mod1, p4.2_mod1, nrow = 2, ncol = 2)

# análisis residuales modelo 2 sin trimestres 2 y 3

p4.1_mod2 <- p4_base + geom_line(aes(y = resid_mod2_sin23trim)) +
  labs(x = 'Tiempo',
       y  = 'Residuales')
p4.1_mod2

p4.2_mod2 <- p4_base + geom_point(aes(x = aju_mod2_sin23trim,
                                      y = resid_mod2_sin23trim)) +
  labs(x = 'Ajustados',
       y = 'Residuales')
p4.2_mod2

p4.3_mod2 <- p4_base + geom_line(aes(y = ivand, colour = 'Real')) + 
  geom_line(aes(y = aju_mod2_sin23trim, colour = 'Ajuste')) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') +
  scale_colour_manual("",
                      values = c('Real'='black', 'Ajuste'='red')) +
  theme(legend.position = c(0.15, 0.85), legend.background = element_blank())
p4.3_mod2

grid.arrange(p4.3_mod2, blankPlot, p4.1_mod2, p4.2_mod2, nrow = 2, ncol = 2)

## 5. Pronósticos para validación cruzada --------------------------------------

# pronósticos modelo 1 sin trimestre 2 y 3
nro_datos_aju <- nrow(datos_ajuste)
t_pron <- (nro_datos_aju + 1):nro_datos
I4_pron <- c(1, 0, 0, 0)
pred_mod1_sin23trim <- predict(mod1_sin23trim,
                               newdata = data.frame(t = t_pron, I4 = I4_pron),
                               interval = 'prediction', level = 0.95)

pred_mod1_sin23trim_esc_orig <-
  exp(pred_mod1_sin23trim) * exp(0.5*sigma_mod1_sin23trim^2)

datos$pron_p_mod1 <- c(rep(NA, nro_datos_aju),
                       pred_mod1_sin23trim_esc_orig[, 1])

datos$pron_ici_mod1 <- c(rep(NA, nro_datos_aju),
                       pred_mod1_sin23trim_esc_orig[, 2])

datos$pron_ics_mod1 <- c(rep(NA, nro_datos_aju),
                         pred_mod1_sin23trim_esc_orig[, 3])

p5_pron_mod1 <- ggplot(data = datos, aes(x = fecha)) + theme_bw() +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = pron_p_mod1, colour = 'Pronóstico'), linetype = 2) +
  geom_line(aes(y = pron_ici_mod1), linetype = 3) +
  geom_line(aes(y = pron_ics_mod1), linetype = 3) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Pronóstico' = 'red')) +
  theme(legend.position = c(0.15, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)
p5_pron_mod1

accuracy(datos$pron_p_mod1[t_pron], datos$ivand[t_pron])

# pronósticos modelo 2 sin trimestre 2 y 3
pred_mod2_sin23trim <- predict(mod2_sin23trim,
                               newdata = data.frame(t = t_pron, I4 = I4_pron),
                               interval = 'prediction', level = 0.95)

datos$pron_p_mod2 <- c(rep(NA, nro_datos_aju), pred_mod2_sin23trim)

p5_pron_mod2 <- ggplot(data = datos, aes(x = fecha)) + theme_bw() +
  geom_line(aes(y = ivand, colour = 'Real')) +
  geom_line(aes(y = pron_p_mod2, colour = 'Pronóstico'), linetype = 2) +
  scale_colour_manual("",
                      values = c('Real' = 'black', 'Pronóstico' = 'red')) +
  theme(legend.position = c(0.15, 0.85), legend.background = element_blank()) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') + ylim(2000, 13000)
p5_pron_mod2

accuracy(datos$pron_p_mod2[t_pron], datos$ivand[t_pron])

grid.arrange(p5_pron_mod1, p5_pron_mod2, nrow = 1, ncol = 2)

## 6. Selección del mejor modelo global, pronóstico ex-ante y estabilidad ------
datos$I2 <- ifelse(datos$tri == 'T2', 1, 0)
datos$I3 <- ifelse(datos$tri == 'T3', 1, 0)
datos$I4 <- ifelse(datos$tri == 'T4', 1, 0)
datos$t <- 1:nro_datos

mod1_full <- lm(log_ivand ~ t + I(t^2) + I4, data = datos)
summary(mod1_full)
qt(0.025, 55) # valor crítico de t

sigma_mod1_full <- summary(mod1_full)$sigma
datos$resid_mod1_full <- resid(mod1_full)
datos$aju_mod1_full <- fitted(mod1_full)
datos$aju_mod1_esc_orig_full <- 
  exp(datos$aju_mod1_full) * exp(0.5*sigma_mod1_full^2)

exp(coef(mod1_full)) * exp(0.5*sigma_mod1_full^2)
resid_mod1_esc_orig_full <- 
  datos$ivand - datos$aju_mod1_esc_orig_full
exp(C_p(resid_mod1_esc_orig_full, type = 'AIC', p = 6))
exp(C_p(resid_mod1_esc_orig_full, type = 'BIC', p = 6))

val_ini_full <- coef(mod1_full)

mod2_full <- nls(ivand ~ exp(b0 + b1*t + b2*I(t^2) + d4*I4),
                      start = list(b0 = val_ini_full[1], b1 = val_ini_full[2],
                                   b2 = val_ini_full[3], d4 = val_ini_full[4]),
                      data = datos)
summary(mod2_full)

resid_mod2_full <- resid(mod2_full)
datos$resid_mod2_full <- resid_mod2_full
datos$aju_mod2_full <- fitted(mod2_full)
exp(C_p(resid_mod2_full, type = 'AIC', p = 6))
exp(C_p(resid_mod2_full, type = 'BIC', p = 6))


p7_base <- ggplot(data = datos, aes(x = fecha)) + theme_bw()

p7.1_mod1 <- p7_base + geom_line(aes(y = resid_mod1_full)) +
  labs(x = 'Tiempo',
       y  = 'Residuales')
p7.1_mod1

p7.2_mod1 <- p7_base + geom_point(aes(x = aju_mod1_full,
                                      y = resid_mod1_full)) +
  labs(x = 'Ajustados',
       y = 'Residuales')
p7.2_mod1

p7.3_mod1 <- p7_base + geom_line(aes(y = ivand, colour = 'Real')) + 
  geom_line(aes(y = aju_mod1_esc_orig_full, colour = 'Ajuste')) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') +
  scale_colour_manual("",
                      values = c('Real'='black', 'Ajuste'='red')) +
  theme(legend.position = c(0.15, 0.85), legend.background = element_blank())
p7.3_mod1

grid.arrange(p7.3_mod1, blankPlot, p7.1_mod1, p7.2_mod1, nrow = 2, ncol = 2)

# análisis residuales modelo 2 sin trimestres 2 y 3

p7.1_mod2 <- p7_base + geom_line(aes(y = resid_mod2_full)) +
  labs(x = 'Tiempo',
       y  = 'Residuales')
p7.1_mod2

p7.2_mod2 <- p7_base + geom_point(aes(x = aju_mod2_full,
                                      y = resid_mod2_full)) +
  labs(x = 'Ajustados',
       y = 'Residuales')
p7.2_mod2

p7.3_mod2 <- p7_base + geom_line(aes(y = ivand, colour = 'Real')) + 
  geom_line(aes(y = aju_mod2_full, colour = 'Ajuste')) +
  labs(x = 'Tiempo',
       y = 'IVA [1000\'s mill de pesos]') +
  scale_colour_manual("",
                      values = c('Real'='black', 'Ajuste'='red')) +
  theme(legend.position = c(0.15, 0.85), legend.background = element_blank())
p7.3_mod2

grid.arrange(p7.3_mod2, blankPlot, p7.1_mod2, p7.2_mod2, nrow = 2, ncol = 2)

pron_mod1_full <- predict(mod1_full,
                          data.frame(t = (nro_datos + 1):(nro_datos + 4),
                                     I4 = c(1, 0, 0, 0)),
                          interval = 'prediction', level = 0.95)
pron_mod1_full_esc_orig <- 
  exp(pron_mod1_full) * exp(0.5*sigma_mod1_full^2)
pron_mod1_full_esc_orig

pron_mod2_full <- predict(mod2_full,
                          data.frame(t = (nro_datos + 1):(nro_datos + 4),
                                     I4 = c(1, 0, 0, 0)),
                          interval = 'prediction', level = 0.95)
pron_mod2_full

S_est <- 4
npar_mod1 <- length(coef(mod1_sin23trim))
k_ini <- max((npar_mod1 + 1), 2 * S_est)
obs_mod_i <- k_ini:nro_datos_aju
nro_mod_ajustados <- length(obs_mod_i)
mtrx_coef <- matrix(data = NA, nrow = nro_mod_ajustados, ncol = npar_mod1*3 + 1)
resid_recur_mod1 <- vector(mode = "numeric", length = nro_mod_ajustados - 1)
for (i in obs_mod_i) {
  mod_tmp <- lm(log_ivand ~ t + I(t^2) + I4,
                data = datos_ajuste[1:i, ])
  sum_tmp <- summary(mod_tmp)
  sigma_mod_tmp <- sum_tmp$sigma
  coef_tmp <- coef(mod_tmp)
  t_val <- qt(0.025, mod_tmp$df.residual)
  lic <- coef_tmp - t_val * sum_tmp$coefficients[, 2]
  lsc <- coef_tmp + t_val * sum_tmp$coefficients[, 2]
  idx_coef <- seq(2, ncol(mtrx_coef), 3)
  idx_lic <- seq(3, ncol(mtrx_coef), 3)
  idx_lsc <- seq(4, ncol(mtrx_coef), 3)
  mtrx_coef[i - k_ini + 1, 1] <- i
  mtrx_coef[i - k_ini + 1, idx_coef] <- coef_tmp
  mtrx_coef[i - k_ini + 1, idx_lic] <- lic
  mtrx_coef[i - k_ini + 1, idx_lsc] <- lsc
  
  if (i < nro_datos_aju) {
    I4_tmp <- 1 * !((i + 1) %% S_est)
    pron_tmp <- predict(mod_tmp,
                        newdata = data.frame(t = i + 1, I4 = I4_tmp),
                        se.fit = TRUE)
    resid_recur_mod1[i - k_ini + 1] <-
      (pron_tmp$fit - datos_ajuste$log_ivand[i + 1]) / (pron_tmp$se.fit*sqrt(i))
  }
}

mtrx_coef
resid_recur_mod1
recresid(mod1_sin23trim, start = 8)

plot(mtrx_coef[,2], type = 'l')
abline(h = coef(mod1_sin23trim)[1], lty = 2)
plot(mtrx_coef[,5], type = 'l')
plot(mtrx_coef[,8], type = 'l')
plot(mtrx_coef[,11], type = 'l')

coef_rec_mod1 <- calc_coef_recursivos(mod1_sin23trim, 8)
coef_rec_mod1
colnames(coef_rec_mod1) <- c('k', 'b0', 'b0.lic', 'b0.lsc', 'b1', 'b1.lic',
                             'b1.lsc', 'b2', 'b2.lic', 'b2.lsc', 'd4',
                             'd4.lic', 'd4.lsc')
coef_rec_mod2 <- calc_coef_recursivos(mod2_sin23trim, 8)
coef_rec_mod2

p7_mod1_base <- ggplot(data = as.data.frame(coef_rec_mod1), aes(x = k)) +
  theme_bw()

p7_mod1_b0 <- p7_mod1_base + geom_line(aes(y = b0)) +
  geom_line(aes(y = b0.lic), linetype = 2) +
  geom_line(aes(y = b0.lsc), linetype = 2) +
  geom_hline(yintercept = coef(mod1_sin23trim)[1], colour = 'red') +
  labs(x = 'n',
       y = expression(hat(beta)[0]))
p7_mod1_b0

p7_mod1_b1 <- p7_mod1_base + geom_line(aes(y = b1)) +
  geom_line(aes(y = b1.lic), linetype = 2) +
  geom_line(aes(y = b1.lsc), linetype = 2) +
  geom_hline(yintercept = coef(mod1_sin23trim)[2], colour = 'red') +
  labs(x = 'n',
       y = expression(hat(beta)[1]))
p7_mod1_b1

p7_mod1_b2 <- p7_mod1_base + geom_line(aes(y = b2)) +
  geom_line(aes(y = b2.lic), linetype = 2) +
  geom_line(aes(y = b2.lsc), linetype = 2) +
  geom_hline(yintercept = coef(mod1_sin23trim)[3], colour = 'red') +
  labs(x = 'n',
       y = expression(hat(beta)[2]))
p7_mod1_b2

p7_mod1_d4 <- p7_mod1_base + geom_line(aes(y = d4)) +
  geom_line(aes(y = d4.lic), linetype = 2) +
  geom_line(aes(y = d4.lsc), linetype = 2) +
  geom_hline(yintercept = coef(mod1_sin23trim)[4], colour = 'red') +
  labs(x = 'n',
       y = expression(hat(delta)[4]))
p7_mod1_d4

grid.arrange(p7_mod1_b0, p7_mod1_b1, p7_mod1_b2, p7_mod1_d4,
             ncol = 2, nrow = 2)

sctest(formula(mod1_sin23trim), data = datos_ajuste, type = 'Rec-CUSUM')
plot(efp(formula(mod1_sin23trim), type="Rec-CUSUM", data = datos_ajuste),
     alpha = 0.05)

p7_mod2_base <- ggplot(data = as.data.frame(coef_rec_mod2), aes(x = k)) +
  theme_bw()

p7_mod2_b0 <- p7_mod2_base + geom_line(aes(y = b0)) +
  geom_line(aes(y = b0.lic), linetype = 2) +
  geom_line(aes(y = b0.lsc), linetype = 2) +
  geom_hline(yintercept = coef(mod2_sin23trim)[1], colour = 'red') +
  labs(x = 'n',
       y = expression(hat(beta)[0]))
p7_mod2_b0

p7_mod2_b1 <- p7_mod2_base + geom_line(aes(y = b1)) +
  geom_line(aes(y = b1.lic), linetype = 2) +
  geom_line(aes(y = b1.lsc), linetype = 2) +
  geom_hline(yintercept = coef(mod2_sin23trim)[2], colour = 'red') +
  labs(x = 'n',
       y = expression(hat(beta)[1]))
p7_mod2_b1

p7_mod2_b2 <- p7_mod2_base + geom_line(aes(y = b2)) +
  geom_line(aes(y = b2.lic), linetype = 2) +
  geom_line(aes(y = b2.lsc), linetype = 2) +
  geom_hline(yintercept = coef(mod2_sin23trim)[3], colour = 'red') +
  labs(x = 'n',
       y = expression(hat(beta)[2]))
p7_mod2_b2

p7_mod2_d4 <- p7_mod2_base + geom_line(aes(y = d4)) +
  geom_line(aes(y = d4.lic), linetype = 2) +
  geom_line(aes(y = d4.lsc), linetype = 2) +
  geom_hline(yintercept = coef(mod2_sin23trim)[4], colour = 'red') +
  labs(x = 'n',
       y = expression(hat(delta)[4]))
p7_mod2_d4

grid.arrange(p7_mod2_b0, p7_mod2_b1, p7_mod2_b2, p7_mod2_d4,
             ncol = 2, nrow = 2)

sctest(formula(mod2_sin23trim), data = datos_ajuste, type = 'Rec-CUSUM')
plot(efp(formula(mod2_sin23trim), type="Rec-CUSUM", data = datos_ajuste), alpha = 0.05)
