#' Cálculo de Tamaño Muestral y Potencia para Estudios de Cohorte
#' 
#' @param alfa Nivel de significación
#' @param potencia Potencia deseada (opcional si se proporciona n_expuestos)
#' @param n_expuestos Número de sujetos expuestos (opcional si se proporciona potencia)
#' @param p0 Incidencia en grupo no expuesto (requerido si metodo = "rho")
#' @param RR Riesgo Relativo a detectar (opcional si se proporciona p1)
#' @param p1 Incidencia en grupo expuesto (opcional si se proporciona RR)
#' @param IC_superior Límite superior del intervalo de confianza del RR
#' @param IC_inferior Límite inferior del intervalo de confianza del RR (opcional si se proporciona EE)
#' @param EE Error estándar del coeficiente (opcional si se proporcionan los IC)
#' @param n_previo Tamaño de muestra del estudio previo (requerido para método "ee")
#' @param r Razón de no expuestos a expuestos (n_no_expuestos / n_expuestos)
#' @param metodo Método de cálculo: "rho" o "ee" (error estándar)
#' @param valores_rho Vector de valores rho para calcular tamaños muestrales (si metodo = "rho")
#' @return Un data frame con resultados según el método elegido
#' @export
SampleRiskRatioMulti <- function(alfa, potencia = NULL, n_expuestos = NULL, 
                                p0 = NULL, RR = NULL, p1 = NULL,
                                IC_superior = NULL, IC_inferior = NULL,
                                EE = NULL, n_previo = NULL, r = 1,
                                metodo = "rho",
                                valores_rho = seq(0, 0.9, by = 0.1)) {
  
  # Verificaciones según el método
  if (metodo == "rho") {
    if (is.null(p0)) {
      stop("Para método 'rho' se requiere p0")
    }
    if (is.null(RR) && is.null(p1)) {
      stop("Para método 'rho' se debe proporcionar RR o p1")
    }
    if (!is.null(RR) && !is.null(p1)) {
      stop("Solo se debe proporcionar RR o p1, no ambos")
    }
  } else if (metodo == "ee") {
    if (is.null(RR)) {
      stop("Para método 'ee' se requiere RR")
    }
    if (is.null(EE) && (is.null(IC_superior) || is.null(IC_inferior))) {
      stop("Para método 'ee' se requiere EE o ambos límites del IC")
    }
    if (is.null(n_previo)) {
      stop("Para método 'ee' se requiere n_previo")
    }
  } else {
    stop("El método debe ser 'rho' o 'ee'")
  }
  
  # Si se proporcionaron los IC, calcular el EE
  if (!is.null(IC_superior) && !is.null(IC_inferior)) {
    EE <- (log(IC_superior) - log(RR)) / qnorm(0.975)
    cat("Error Estándar calculado:", EE, "\n")
  }
  
  if (metodo == "rho") {
    # Calcular p1 si se proporcionó RR
    if (!is.null(RR)) {
      p1 <- p0 * RR
    } else {
      RR <- p1/p0
    }
    
    z_alfa <- qnorm(1 - alfa/2)
    
    if (!is.null(potencia)) {
      # Calcular tamaño muestral
      z_beta <- qnorm(potencia)
      
      calcular_n <- function(rho) {
        numerador <- (z_alfa * sqrt((1 + 1/r) * p0 * (1 - p0)) + 
                       z_beta * sqrt(p1 * (1 - p1) + p0 * (1 - p0) / r))^2
        denominador <- (p1 - p0)^2 * (1 - rho^2)
        
        n_expuestos <- ceiling(numerador / denominador)
        n_no_expuestos <- ceiling(n_expuestos * r)
        
        return(c(n_expuestos, n_no_expuestos))
      }
      
      resultados <- data.frame(
        rho = valores_rho,
        n_expuestos = sapply(valores_rho, function(rho) calcular_n(rho)[1]),
        n_no_expuestos = sapply(valores_rho, function(rho) calcular_n(rho)[2])
      )
      
      resultados$tamaño_muestral_total <- resultados$n_expuestos + resultados$n_no_expuestos
    } else {
      # Calcular potencia
      n_no_expuestos <- ceiling(n_expuestos * r)
      
      calcular_potencia <- function(rho) {
        numerador <- sqrt(n_expuestos * (p1 - p0)^2 * (1 - rho^2)) - 
          z_alfa * sqrt((1 + 1/r) * p0 * (1 - p0))
        denominador <- sqrt(p1 * (1 - p1) + p0 * (1 - p0) / r)
        
        potencia <- pnorm(numerador / denominador)
        return(potencia)
      }
      
      resultados <- data.frame(
        rho = valores_rho,
        potencia = sapply(valores_rho, calcular_potencia)
      )
      
      resultados$n_expuestos <- n_expuestos
      resultados$n_no_expuestos <- n_no_expuestos
      resultados$tamaño_muestral_total <- n_expuestos + n_no_expuestos
    }
  } else {
    # Método basado en error estándar
    z_alfa <- qnorm(1 - alfa/2)
    if (!is.null(potencia)) {
      z_gamma <- qnorm(potencia)
      n <- ceiling((z_alfa + z_gamma)^2 * n_previo * EE^2 / (log(RR))^2)
      resultados <- data.frame(
        tamaño_muestral_total = n
      )
      resultados$n_expuestos <- ceiling(n / (1 + r))
      resultados$n_no_expuestos <- ceiling(n * r / (1 + r))
    } else {
      z_gamma <- sqrt(n_expuestos * (log(RR))^2 / (n_previo * EE^2)) - z_alfa
      potencia <- pnorm(z_gamma)
      resultados <- data.frame(
        potencia = potencia
      )
      resultados$n_expuestos <- n_expuestos
      resultados$n_no_expuestos <- ceiling(n_expuestos * r)
      resultados$tamaño_muestral_total <- n_expuestos + resultados$n_no_expuestos
    }
  }
  return(resultados)
}

#' Cálculo Logístico para Estudios de Cohorte
#'
#' @param n_expuestos Número de sujetos expuestos requeridos
#' @param n_no_expuestos Número de sujetos no expuestos requeridos
#' @param tasa_reclutamiento Tasa de reclutamiento por día
#' @param tasa_rechazo Tasa de rechazo esperada
#' @param tasa_perdida_seguimiento Tasa esperada de pérdida en el seguimiento
#' @param tasa_elegibilidad Tasa de elegibilidad esperada
#' @param dias_laborables_mes Número de días laborables por mes
#' @return Una lista con cálculos logísticos para el estudio
#' @export
logistica_estudio_cohorte <- function(n_expuestos, 
                                     n_no_expuestos,
                                     tasa_reclutamiento,
                                     tasa_rechazo,
                                     tasa_perdida_seguimiento,
                                     tasa_elegibilidad,
                                     dias_laborables_mes) {
  
  tamaño_muestral_total <- n_expuestos + n_no_expuestos
  
  # Ajuste por pérdidas en el seguimiento
  muestra_con_perdidas <- tamaño_muestral_total / (1 - tasa_perdida_seguimiento)
  
  # Ajuste por rechazos
  muestra_con_rechazos <- muestra_con_perdidas / (1 - tasa_rechazo)
  
  # Ajuste por elegibilidad
  muestra_a_evaluar <- muestra_con_rechazos / tasa_elegibilidad
  
  # Cálculo de tiempos
  dias_reclutamiento <- ceiling(muestra_a_evaluar / tasa_reclutamiento)
  meses_reclutamiento <- dias_reclutamiento / dias_laborables_mes
  
  # Resultados
  resultados <- list(
    muestra_final = tamaño_muestral_total,
    muestra_con_perdidas = ceiling(muestra_con_perdidas),
    muestra_para_enrolar = ceiling(muestra_con_rechazos),
    muestra_a_evaluar = ceiling(muestra_a_evaluar),
    dias_reclutamiento = dias_reclutamiento,
    meses_reclutamiento = round(meses_reclutamiento, 2)
  )
  
  # Imprimir resumen
  cat("\nResumen de logística del estudio:\n")
  cat("----------------------------------------\n")
  cat("Muestra final requerida:", tamaño_muestral_total, "\n")
  cat("Muestra considerando pérdidas:", ceiling(muestra_con_perdidas),
      "(", tasa_perdida_seguimiento*100, "% de pérdidas)\n")
  cat("Muestra a enrolar:", ceiling(muestra_con_rechazos),
      "(", tasa_rechazo*100, "% de rechazo)\n")
  cat("Muestra a evaluar:", ceiling(muestra_a_evaluar),
      "(", tasa_elegibilidad*100, "% de elegibilidad)\n")
  cat("Días necesarios:", dias_reclutamiento,
      "(", tasa_reclutamiento, "personas por día)\n")
  cat("Meses necesarios:", round(meses_reclutamiento, 2),
      "(", dias_laborables_mes, "días laborables por mes)\n")
  cat("----------------------------------------\n")
  
  return(invisible(resultados))
}
