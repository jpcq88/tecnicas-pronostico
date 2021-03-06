#La siguiente linea permite leer los datos guardados en archivo
#de texto en disco local. file.choose() abre ventana de windows para
#explorar y ubicar el archivo. header=T para indicar que los datos
#est�n encabezados por nombre de la columna; si no hay encabezado, se toma header=F
datos.indpro=scan(file.choose(),dec = ".")
datos.indpro=ts(datos.indpro,frequency=12,start=c(1950,1))
datos.indpro
win.graph()
plot(datos.indpro)

plot(decompose(datos.indpro,type="additive"))

win.graph(width=5,height=4,pointsize=8)
boxplot(datos.indpro~cycle(datos.indpro),names=month.abb)

#Creando �ndice t de tiempo
t=1:length(datos.indpro)

#Ajuste modelo de tendencia cuadr�tica
modelo1=lm(datos.indpro~t+I(t^2))
summary(modelo1)
ajuste=ts(fitted(modelo1),frequency=12,start=c(1950,1))

#Gr�fico de la serie y su ajuste
win.graph()
plot(datos.indpro)
lines(ajuste,col=2)
legend("topleft",legend=c("Serie real","Serie ajustada modelo 1"),col=c(1,2),lty=1)


# Gr�fico de residuales comunes  versus variable respuesta
win.graph()
plot(fitted(modelo1), residuals(modelo1), type="p", main="Ajustados vs.Residules comunes modelo 1")
abline(h=0,col=2)

# Gr�fico de residuales comunes  versus tiempo
win.graph()
plot(t, residuals(modelo1), type="l", main="Tiempo vs. Residuales comunes modelo 1")
abline(h=0,col=2)

# Prueba de Normalidad
shapiro.test(residuals(modelo1))
win.graph()
qqnorm(residuals(modelo1),main="Gr�fico de normalidad en modelo 1")
qqline(residuals(modelo1),col=2)



#Ajuste modelo de tendencia c�bica
modelo2=lm(datos.indpro~t+I(t^2)+I(t^3))
summary(modelo2)
ajuste2=ts(fitted(modelo2),frequency=12,start=c(1950,1))

#Gr�fico de la serie y su ajuste
win.graph()
plot(datos.indpro)
lines(ajuste2,col=2)
legend("topleft",legend=c("Serie real","Serie ajustada modelo 2"),col=c(1,2),lty=1)

# Gr�fico de residuales comunes  versus variable respuesta
win.graph()
plot(fitted(modelo2), residuals(modelo2), type="p", main="Ajustados vs.Residules comunes modelo 2")
abline(h=0,col=2)

# Gr�fico de residuales comunes  versus tiempo
win.graph()
plot(t, residuals(modelo2), type="l", main="Tiempo vs. Residuales comunes modelo 2")
abline(h=0,col=2)

# Prueba de Normalidad
shapiro.test(residuals(modelo2))
win.graph()
qqnorm(residuals(modelo2), main="Gr�fico de normalidad en modelo 2")
qqline(residuals(modelo2),col=2)


#Ajuste modelo auxiliar log-c�bico para hallar valores iniciales
#de par�metros en modelo exponencia-c�bico
aux=lm(log(datos.indpro)~t+I(t^2)+I(t^3))
betas.ini=coef(aux) #guardando par�metros estimados modelo auxiliar

#Ajuste modelo exponencial-c�bico
modelo3=nls(datos.indpro~exp(b0+b1*t+b2*I(t^2)+b3*I(t^3)),start=list(b0=betas.ini[1],b1=betas.ini[2],b2=betas.ini[3],b3=betas.ini[4]))
ajuste3=ts(fitted(modelo3),frequency=12,start=c(1950,1))

#Gr�fico de la serie y su ajuste
win.graph()
plot(datos.indpro)
lines(ajuste3,col=2)
legend("topleft",legend=c("Serie real","Serie ajustada modelo 3"),col=c(1,2),lty=1)

# Gr�fico de residuales comunes  versus variable respuesta
win.graph()
plot(fitted(modelo3), residuals(modelo3), type="p", main="Ajustados vs.Residules comunes modelo 3")
abline(h=0,col=2)

# Gr�fico de residuales comunes  versus tiempo
win.graph()
plot(t, residuals(modelo3), type="l", main="Tiempo vs. Residuales comunes modelo 3")
abline(h=0,col=2)

# Prueba de Normalidad
shapiro.test(residuals(modelo3))
win.graph()
qqnorm(residuals(modelo3), main="Gr�fico de normalidad en modelo 3")
qqline(residuals(modelo3),col=2)

AIC(modelo1)
AIC(modelo2)
AIC(modelo3)
BIC(modelo1)
BIC(modelo2)
BIC(modelo3)

#Ajuste de modelos excluyendo los tres primeros a�osrm(list=ls(all=TRUE)) #elimina todos los objetos creados
datos.indpro=scan(file.choose(),dec = ".") #lea de nuevo los datos

#Datos desde enero de 1953
datos.indpro2=ts(datos.indpro[37:length(datos.indpro)],freq=12,start=c(1953,1))

win.graph()
plot(datos.indpro2)

#Creando �ndice t de tiempo
t=1:length(datos.indpro2)

#Ajuste modelo de tendencia cuadr�tica
modelo1b=lm(datos.indpro2~t+I(t^2))
summary(modelo1b)
ajuste1b=ts(fitted(modelo1b),frequency=12,start=c(1953,1))

#Gr�fico de la serie y su ajuste
win.graph()
plot(datos.indpro2)
lines(ajuste1b,col=2)
legend("topleft",legend=c("Serie real","Serie Ajustada modelo 1"),col=c(1,2),lty=1)

# Gr�fico de residuales comunes  versus variable respuesta
win.graph()
plot(fitted(modelo1b), residuals(modelo1b), type="p", main="Ajustados vs.Residules comunes modelo 1\nDatos desde 1953.I")
abline(h=0,col=2)

# Gr�fico de residuales comunes  versus tiempo
win.graph()
plot(t, residuals(modelo1b), type="l", main="Tiempo vs. Residuales comunes modelo 1\nDatos desde 1953.I")
abline(h=0,col=2)

# Prueba de Normalidad
shapiro.test(residuals(modelo1b))
win.graph()
qqnorm(residuals(modelo1b),main="Gr�fico de normalidad en modelo 1\nDatos desde 1953.I")
qqline(residuals(modelo1b),col=2)

#Correr modelo c�bico y el exponencial-c�bico para estos datos y comparar con el modelo cuadr�tico con AIC, BIC y comportamiento de residuos
AIC(modelo1b)
BIC(modelo1b)
