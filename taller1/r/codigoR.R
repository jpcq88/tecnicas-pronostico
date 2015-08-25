#------------------------------------------------------------------------------#
# Código R taller 1: Análisis índice de productividad Canadá                   #
#------------------------------------------------------------------------------#

# Librerías requeridas
library(ggplot2) # graficos
library(grid) # para poder utilizar la función unit() en ggplot2
library(gridExtra) # funciones auxiliares graficar grid.arrange()
library(dplyr) # manipulación de datos
library(tidyr) # organizar los datos
library(xtable) # exportar tablas en formato LaTeX
theme_set(theme_bw(base_size = 7)) # definir tamaño de fuente gráficos ggplot2

## Funciones auxiliares
ggqqplot <- function(mod){
  # Devuelve un gráfico qqplot al estilo ggplot2
  
  # res_std <- rstudent(mod)
  res_std <- resid(mod) # descomentar para usar con modelo no lineal
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

## Lectura de datos
datos_orig <- read.table(file = '../datos/IndProductividad.txt', header = FALSE)
head(datos_orig)

## Organizamos los datos y le damos el formato adecuado
ind_prod <- gather(as.data.frame(t(datos_orig)))[, 2]
fecha <- seq.Date(from = as.Date('1950-01-01', '%Y-%m-%d'),
                  to = as.Date('1973-12-01', '%Y-%m-%d'),
                  by = 'month')
mes <- reorder(months.Date(fecha), rep(1:12, 24))
datos <- data.frame(mes, fecha, ind_prod)
head(datos)
str(datos)

## Análisis exploratorio inicial
summary(datos)

## Punto 1
p1_ind_prod <- ggplot(data = datos) + theme_bw() +
  geom_line(aes(x = fecha, y = ind_prod)) +
  labs(x = '', y = 'Índice de pruductividad') +
  theme(
    axis.line = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )
p1_ind_prod # ver el gráfico

# Guardar gráfico
ggsave(filename = '../img/p1_ind_prod.pdf', plot = p1_ind_prod,
       width = 7, height = 2.5)

# Punto 1
datos_ts <- ts(datos$ind_prod, start = c(1950, 1), end = c(1973, 12),
               frequency = 12)
ts_desc <- decompose(datos_ts)
plot(ts_desc)
str(ts_desc)

pdf('../img/p1_descomp.pdf', width = 7, height = 4)
opar <- par()
par(mfrow = c(3, 1),
    oma = c(1, 1, 1, 1),
    mar = c(1, 4, 1, 1),
    xaxt = 'n',
    cex.axis = 0.6,
    cex.lab = 0.9,
    lwd = 0.5)
plot(ts_desc$trend,
     xlab = '',
     ylab = expression(paste(T[t])))
plot(ts_desc$seasonal,
     xlab = '',
     ylab = expression(paste(S[t])))
abline(v = 1950:1974, lty = 2, col =  'red')
par(xaxt = 's')
plot(ts_desc$random,
     xlab = 'Tiempo',
     ylab = expression(paste(E[t])))
abline(h = 0, lty = 2, col =  'red')
par(opar)
dev.off()

prom_mes <- summarise(group_by(datos, mes), prom = mean(ind_prod))

meses <- c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun',
           'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic')

p1_boxplot_mes <- ggplot(data = datos) + theme_bw() +
  geom_boxplot(aes(x = mes, y = ind_prod)) +
  geom_line(data = prom_mes, aes(x = 1:12, y = prom),
            size = 0.2, colour = 'red', linetype = 6) +
  scale_x_discrete(labels = meses) +
  labs(x = '', y = 'Índice de pruductividad') +
  theme(
    axis.line = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

p1_boxplot_mes
ggsave(filename = '../img/p1_boxplot_mes.pdf', plot = p1_boxplot_mes,
       width = 4.5, height = 3)

## Punto 2

tiempo <- 1:length(datos_ts)
mod1 <- lm(datos_ts ~ tiempo + I(tiempo^2))
summary(mod1)
xtable(summary(mod1),
       caption = 'Ajuste modelo tendencia cuadrática',
       label = 'tab:mod1_comp',
       align = 'lrrrr')

datos$aju_mod1_com <- fitted(mod1)
datos$res_mod1_com <- rstudent(mod1)

p2_mod1_g1 <- ggplot(data = datos) + theme_bw(8) +
  geom_line(aes(x = fecha, y = ind_prod)) +
  geom_line(aes(x = fecha, y = aju_mod1_com), colour = 'red') +
  labs(x = '',
       y = 'Índice de producción',
       title = 'Serie original vs modelo cuadrático ajustado')

p2_mod1_g2 <- ggplot(data = datos) + theme_bw(8) +
  geom_point(aes(x = aju_mod1_com, y = res_mod1_com), size = 1) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = 'Valores ajustados',
       y = 'Residuales estudentizados',
       title = 'Análisis residuales estudentizados')

p2_mod1_g3 <- ggplot(data = datos) + theme_bw(8) +
  geom_line(aes(x = fecha, y = res_mod1_com)) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = '',
       y = 'Residuales estudentizados',
       title = 'Análisis residuales estudentizados')

p2_mod1_g4 <- ggqqplot(mod1)

pdf('../img/p2_diag_mod1_comp.pdf', width = 7, height = 7)
grid.arrange(p2_mod1_g1, p2_mod1_g2, p2_mod1_g3, p2_mod1_g4,
             ncol = 2, nrow = 2)
dev.off()

mod2 <- lm(datos_ts ~ tiempo + I(tiempo^2) + I(tiempo^3))
summary(mod2)
xtable(summary(mod2),
       caption = 'Ajuste modelo tendencia cúbica',
       label = 'tab:mod2_comp',
       align = 'lrrrr')

datos$aju_mod2_com <- fitted(mod2)
datos$res_mod2_com <- rstudent(mod2)

p2_mod2_g1 <- ggplot(data = datos) + theme_bw(8) +
  geom_line(aes(x = fecha, y = ind_prod)) +
  geom_line(aes(x = fecha, y = aju_mod2_com), colour = 'red') +
  labs(x = '',
       y = 'Índice de producción',
       title = 'Serie original vs modelo cúbico ajustado')

p2_mod2_g2 <- ggplot(data = datos) + theme_bw(8) +
  geom_point(aes(x = aju_mod2_com, y = res_mod2_com), size = 1) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = 'Valores ajustados',
       y = 'Residuales estudentizados',
       title = 'Análisis residuales estudentizados')

p2_mod2_g3 <- ggplot(data = datos) + theme_bw(8) +
  geom_line(aes(x = fecha, y = res_mod2_com)) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = '',
       y = 'Residuales estudentizados',
       title = 'Análisis residuales estudentizados')

p2_mod2_g4 <- ggqqplot(mod2)

pdf('../img/p2_diag_mod2_comp.pdf', width = 7, height = 7)
grid.arrange(p2_mod2_g1, p2_mod2_g2, p2_mod2_g3, p2_mod2_g4,
             ncol = 2, nrow = 2)
dev.off()

mod3_init <- lm(log(datos_ts) ~ tiempo + I(tiempo^2) + I(tiempo^3))
coef_init <- mod3_init$coefficients

mod3 <- nls(datos_ts ~ exp(b0 + b1*tiempo + b2*I(tiempo^2) + b3*I(tiempo^3)),
            start = list(b0 = coef_init[1],
                         b1 = coef_init[2],
                         b2 = coef_init[3],
                         b3 = coef_init[4]))
summary(mod3)

datos$aju_mod3_com <- fitted(mod3)
datos$res_mod3_com <- resid(mod3)

p2_mod3_g1 <- ggplot(data = datos) + theme_bw(8) +
  geom_line(aes(x = fecha, y = ind_prod)) +
  geom_line(aes(x = fecha, y = aju_mod3_com), colour = 'red') +
  labs(x = '',
       y = 'Índice de producción',
       title = 'Serie original vs modelo cúbico expoencial ajustado')

p2_mod3_g2 <- ggplot(data = datos) + theme_bw(8) +
  geom_point(aes(x = aju_mod3_com, y = res_mod3_com), size = 1) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = 'Valores ajustados',
       y = 'Residuales',
       title = 'Análisis residuales')

p2_mod3_g3 <- ggplot(data = datos) + theme_bw(8) +
  geom_line(aes(x = fecha, y = res_mod3_com)) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = '',
       y = 'Residuales',
       title = 'Análisis residuales')

p2_mod3_g4 <- ggqqplot(mod3)

pdf('../img/p2_diag_mod3_comp.pdf', width = 7, height = 7)
grid.arrange(p2_mod3_g1, p2_mod3_g2, p2_mod3_g3, p2_mod3_g4,
             ncol = 2, nrow = 2)
dev.off()

# Comparación de los modelos
summary(mod1)
summary(mod2)
summary(mod3)
AIC(mod1)
AIC(mod2)
AIC(mod3)
BIC(mod1)
BIC(mod2)
BIC(mod3)

## Punto 3

tiempo2 <- 1:(length(datos_ts)-3*12)
datos_ts2 <- datos_ts[37:288]
mod12 <- lm(datos_ts2 ~ tiempo2 + I(tiempo2^2))
summary(mod1)
xtable(summary(mod1),
       caption = 'Ajuste modelo tendencia cuadrática sin tres primeros años',
       label = 'tab:mod1_sin3',
       align = 'lrrrr')

datos2 <- datos[37:288,]
datos2$aju_mod12_com <- fitted(mod12)
datos2$res_mod12_com <- rstudent(mod12)

p3_mod12_g1 <- ggplot(data = datos2) + theme_bw(8) +
  geom_line(aes(x = fecha, y = ind_prod)) +
  geom_line(aes(x = fecha, y = aju_mod12_com), colour = 'red') +
  labs(x = '',
       y = 'Índice de producción',
       title = 'Serie original vs modelo cuadrático ajustado')

p3_mod12_g2 <- ggplot(data = datos2) + theme_bw(8) +
  geom_point(aes(x = aju_mod12_com, y = res_mod12_com), size = 1) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = 'Valores ajustados',
       y = 'Residuales estudentizados',
       title = 'Análisis residuales estudentizados')

p3_mod12_g3 <- ggplot(data = datos2) + theme_bw(8) +
  geom_line(aes(x = fecha, y = res_mod12_com)) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = '',
       y = 'Residuales estudentizados',
       title = 'Análisis residuales estudentizados')

p3_mod12_g4 <- ggqqplot(mod12)

pdf('../img/p3_diag_mod1_sin3.pdf', width = 7, height = 7)
grid.arrange(p3_mod12_g1, p3_mod12_g2, p3_mod12_g3, p3_mod12_g4,
             ncol = 2, nrow = 2)
dev.off()

mod22 <- lm(datos_ts2 ~ tiempo2 + I(tiempo2^2) + I(tiempo2^3))
summary(mod22)
xtable(summary(mod22),
       caption = 'Ajuste modelo tendencia cúbica sin tres primeros años',
       label = 'tab:mod2_sin3',
       align = 'lrrrr')

datos2$aju_mod22_com <- fitted(mod22)
datos2$res_mod22_com <- rstudent(mod22)

p3_mod22_g1 <- ggplot(data = datos2) + theme_bw(8) +
  geom_line(aes(x = fecha, y = ind_prod)) +
  geom_line(aes(x = fecha, y = aju_mod22_com), colour = 'red') +
  labs(x = '',
       y = 'Índice de producción',
       title = 'Serie original vs modelo cúbico ajustado')

p3_mod22_g2 <- ggplot(data = datos2) + theme_bw(8) +
  geom_point(aes(x = aju_mod22_com, y = res_mod22_com), size = 1) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = 'Valores ajustados',
       y = 'Residuales estudentizados',
       title = 'Análisis residuales estudentizados')

p3_mod22_g3 <- ggplot(data = datos2) + theme_bw(8) +
  geom_line(aes(x = fecha, y = res_mod22_com)) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = '',
       y = 'Residuales estudentizados',
       title = 'Análisis residuales estudentizados')

p3_mod22_g4 <- ggqqplot(mod22)

pdf('../img/p3_diag_mod2_sin3.pdf', width = 7, height = 7)
grid.arrange(p3_mod22_g1, p3_mod22_g2, p3_mod22_g3, p3_mod22_g4,
             ncol = 2, nrow = 2)
dev.off()

mod32_init <- lm(log(datos_ts2) ~ tiempo2 + I(tiempo2^2) + I(tiempo2^3))
coef_init2 <- mod32_init$coefficients

mod32 <- nls(datos_ts2 ~ exp(b0+b1*tiempo2 + b2*I(tiempo2^2) + b3*I(tiempo2^3)),
            start = list(b0 = coef_init2[1],
                         b1 = coef_init2[2],
                         b2 = coef_init2[3],
                         b3 = coef_init2[4]))
summary(mod32)

datos2$aju_mod32_com <- fitted(mod32)
datos2$res_mod32_com <- resid(mod32)

p3_mod32_g1 <- ggplot(data = datos2) + theme_bw(8) +
  geom_line(aes(x = fecha, y = ind_prod)) +
  geom_line(aes(x = fecha, y = aju_mod32_com), colour = 'red') +
  labs(x = '',
       y = 'Índice de producción',
       title = 'Serie original vs modelo cúbico expoencial ajustado')

p3_mod32_g2 <- ggplot(data = datos2) + theme_bw(8) +
  geom_point(aes(x = aju_mod32_com, y = res_mod32_com), size = 1) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = 'Valores ajustados',
       y = 'Residuales',
       title = 'Análisis residuales')

p3_mod32_g3 <- ggplot(data = datos2) + theme_bw(8) +
  geom_line(aes(x = fecha, y = res_mod32_com)) +
  geom_abline(intercept = 0, slope = 0, colour = 'red', linetype = 2) +
  labs(x = '',
       y = 'Residuales',
       title = 'Análisis residuales')

p3_mod32_g4 <- ggqqplot(mod32)

pdf('../img/p3_diag_mod3_sin3.pdf', width = 7, height = 7)
grid.arrange(p3_mod32_g1, p3_mod32_g2, p3_mod32_g3, p3_mod32_g4,
             ncol = 2, nrow = 2)
dev.off()

# Comparación de los modelos
summary(mod12)
summary(mod22)
summary(mod32)
AIC(mod12)
AIC(mod22)
AIC(mod32)
BIC(mod12)
BIC(mod22)
BIC(mod32)
