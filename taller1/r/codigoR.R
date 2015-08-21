#------------------------------------------------------------------------------#
# Código R taller 1: Análisis índice de productividad Canadá                   #
#------------------------------------------------------------------------------#

# Librerías requeridas
library(ggplot2) # graficos
theme_set(theme_bw(base_size = 7)) # definir tamaño de fuente gráficos
library(dplyr) # manipulación de datos
library(tidyr) # organizar los datos

# Lectura de datos
datos_orig <- read.table(file = '../datos/IndProductividad.txt', header = FALSE)
head(datos_orig)

# Organizamos los datos y le damos el formato adecuado
ind_prod <- gather(as.data.frame(t(datos_orig)))[, 2]
fecha <- seq.Date(from = as.Date('1950-01-01', '%Y-%m-%d'),
                  to = as.Date('1973-12-01', '%Y-%m-%d'),
                  by = 'month')
mes <- reorder(months.Date(fecha), rep(1:12, 24))
datos <- data.frame(mes, fecha, ind_prod)
head(datos)
str(datos)

# Análisis exploratorio inicial
summary(datos)

# crear gráfico exploratorio y almacenarlo en p1_ind_prod
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
ggsave(filename = '../img/p1_ind_prod.pdf', p1_ind_prod, width = 7, height = 4)

# Punto 1
datos_ts <- ts(datos$ind_prod, start = c(1950, 1), end = c(1973, 12),
               frequency = 12)
ts_desc <- decompose(datos_ts)
plot(ts_desc)
str(ts_desc)

opar <- par()
par(mfrow = c(4, 1))
plot(ts_desc$x, mar = c(0.1, 0.1, 0.1, 0.1), oma = c(0.1, 0.1, 0.1, 0.1),
     xlab = '')
par(opar)

summarise(group_by(datos, mes), mean(ind_prod))

meses <- c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun',
           'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic')

p1_boxplot_mes <- ggplot(data = datos) + theme_bw() +
  geom_boxplot(aes(x = mes, y = ind_prod)) +
  scale_x_discrete(labels = meses) +
  labs(x = '', y = 'Índice de pruductividad') +
  theme(
    axis.line = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

p1_boxplot_mes
ggsave(filename = '../img/p1_boxplot_mes.pdf', p1_boxplot_mes,
       width = 7, height = 4)
