calcular_entropia_janela_deslizante_ecg <- function(ecg_df, janela_tamanho = 30) {
  if (!exists("Entropy.P")) {
    source("EntropyEstimationParetoKernel2.R")
  }
  
  if (!all(c("time", "mV") %in% names(ecg_df))) {
    stop("O dataframe precisa ter as colunas 'time' e 'mV'.")
  }

  if (janela_tamanho <= 0) {
    stop("O tamanho da janela precisa ser maior que zero.")
  }

  if (nrow(ecg_df) < janela_tamanho) {
    stop("O dataframe tem menos linhas que o tamanho da janela.")
  }
  
  inicio_blocos <- seq(1, nrow(ecg_df) - janela_tamanho + 1, by = janela_tamanho)
  fim_blocos <- inicio_blocos + janela_tamanho - 1
  tamanho_intervalos <- length(inicio_blocos)
  
  blocos_ecg <- data.frame(
    bloco = seq_len(tamanho_intervalos),
    inicio = ecg_df$time[inicio_blocos],
    fim = ecg_df$time[fim_blocos]
  )
  
  sim_ecg_heavy <- numeric(tamanho_intervalos)
  alpha_ecg_heavy <- numeric(tamanho_intervalos)
  
  sim_ecg_light <- numeric(tamanho_intervalos)
  alpha_ecg_light <- numeric(tamanho_intervalos)
  
  sim_ecg_dirac <- numeric(tamanho_intervalos)
  alpha_ecg_dirac <- numeric(tamanho_intervalos)
  
  tipo_ecg_melhor <- character(tamanho_intervalos)
  
  for (b in seq_len(tamanho_intervalos)) {
    
    vetor_ecg <- ecg_df$mV[inicio_blocos[b]:fim_blocos[b]]
    
    h.P.sim_heavy <- optim(
      par = c(0.5, 1),
      fn = function(par) Entropy.P(x.var = vetor_ecg, par = par)
    )
    
    h.P.sim_light <- optim(
      par = c(5, 1),
      fn = function(par) Entropy.P(x.var = vetor_ecg, par = par)
    )
    
    h.P.sim_dirac <- optim(
      par = c(1000, 1),
      fn = function(par) Entropy.P(x.var = vetor_ecg, par = par)
    )
    
    sim_ecg_heavy[b] <- h.P.sim_heavy$value
    alpha_ecg_heavy[b] <- h.P.sim_heavy$par[1]
    
    sim_ecg_light[b] <- h.P.sim_light$value
    alpha_ecg_light[b] <- h.P.sim_light$par[1]
    
    sim_ecg_dirac[b] <- h.P.sim_dirac$value
    alpha_ecg_dirac[b] <- h.P.sim_dirac$par[1]
    
    valores <- c(
      heavy = h.P.sim_heavy$value,
      light = h.P.sim_light$value,
      dirac = h.P.sim_dirac$value
    )
    
    tipo_ecg_melhor[b] <- names(valores)[which.min(valores)]
  }
  
  blocos_ecg$sim_heavy <- sim_ecg_heavy
  blocos_ecg$alpha_heavy <- alpha_ecg_heavy
  blocos_ecg$sim_light <- sim_ecg_light
  blocos_ecg$alpha_light <- alpha_ecg_light
  blocos_ecg$sim_dirac <- sim_ecg_dirac
  blocos_ecg$alpha_dirac <- alpha_ecg_dirac
  blocos_ecg$tipo_melhor <- tipo_ecg_melhor
  
  blocos_ecg
}


calcular_entropia_janela_deslizante_ecg_alpha_otimo <- function(
    ecg_df,
    janela_tamanho = 30,
    par_inicial = c(alpha = 5, beta = 1),
    alpha_dirac_limite = 1000
) {
  if (!exists("Entropy.P")) {
    source("EntropyEstimationParetoKernel2.R")
  }
  
  if (!all(c("time", "mV") %in% names(ecg_df))) {
    stop("O dataframe precisa ter as colunas 'time' e 'mV'.")
  }

  if (janela_tamanho <= 0) {
    stop("O tamanho da janela precisa ser maior que zero.")
  }

  if (nrow(ecg_df) < janela_tamanho) {
    stop("O dataframe tem menos linhas que o tamanho da janela.")
  }
  
  inicio_blocos <- seq(1, nrow(ecg_df) - janela_tamanho + 1, by = janela_tamanho)
  fim_blocos <- inicio_blocos + janela_tamanho - 1
  tamanho_intervalos <- length(inicio_blocos)
  
  blocos_ecg <- data.frame(
    bloco = seq_len(tamanho_intervalos),
    inicio = ecg_df$time[inicio_blocos],
    fim = ecg_df$time[fim_blocos]
  )
  
  sim_ecg_otimo <- numeric(tamanho_intervalos)
  alpha_ecg_otimo <- numeric(tamanho_intervalos)
  beta_ecg_otimo <- numeric(tamanho_intervalos)
  tipo_ecg_melhor <- character(tamanho_intervalos)
  convergencia <- integer(tamanho_intervalos)
  
  classificar_alpha <- function(alpha) {
    if (is.infinite(alpha) || alpha >= alpha_dirac_limite) {
      return("dirac")
    }
    
    if (alpha <= 2) {
      return("heavy")
    }
    
    "light"
  }
  
  for (b in seq_len(tamanho_intervalos)) {
    vetor_ecg <- ecg_df$mV[inicio_blocos[b]:fim_blocos[b]]
    penalidade <- 1e100
    
    objetivo <- function(par) {
      alpha <- par[1]
      beta <- par[2]
      
      if (!is.finite(alpha) || !is.finite(beta) || alpha <= 0 || beta <= 0) {
        return(penalidade)
      }
      
      valor <- tryCatch(
        Entropy.P(x.var = vetor_ecg, par = c(alpha, beta)),
        error = function(e) penalidade
      )
      
      if (!is.finite(valor)) {
        return(penalidade)
      }
      
      valor
    }
    
    h.P.sim_otimo <- optim(
      par = par_inicial,
      fn = objetivo,
      method = "L-BFGS-B",
      lower = c(1e-8, 1e-8),
      upper = c(alpha_dirac_limite, Inf)
    )
    
    sim_ecg_otimo[b] <- h.P.sim_otimo$value
    alpha_ecg_otimo[b] <- h.P.sim_otimo$par[1]
    beta_ecg_otimo[b] <- h.P.sim_otimo$par[2]
    tipo_ecg_melhor[b] <- classificar_alpha(h.P.sim_otimo$par[1])
    convergencia[b] <- h.P.sim_otimo$convergence
  }
  
  blocos_ecg$sim_otimo <- sim_ecg_otimo
  blocos_ecg$alpha_otimo <- alpha_ecg_otimo
  blocos_ecg$beta_otimo <- beta_ecg_otimo
  blocos_ecg$tipo_melhor <- tipo_ecg_melhor
  blocos_ecg$convergencia <- convergencia
  
  blocos_ecg
}


calcular_entropia_janela_deslizante_ecg_com_sobreposicao <- function(ecg_df, janela_tamanho = 30) {
  if (!exists("Entropy.P")) {
    source("EntropyEstimationParetoKernel2.R")
  }

  if (!all(c("time", "mV") %in% names(ecg_df))) {
    stop("O dataframe precisa ter as colunas 'time' e 'mV'.")
  }

  if (janela_tamanho <= 0) {
    stop("O tamanho da janela precisa ser maior que zero.")
  }

  tamanho_intervalos <- nrow(ecg_df) - janela_tamanho + 1

  if (tamanho_intervalos < 1) {
    stop("O dataframe tem menos linhas que o tamanho da janela.")
  }

  blocos_ecg <- data.frame(
    bloco = seq_len(tamanho_intervalos),
    inicio = ecg_df$time[seq_len(tamanho_intervalos)],
    fim = ecg_df$time[seq_len(tamanho_intervalos) + janela_tamanho - 1]
  )

  sim_ecg_heavy <- numeric(tamanho_intervalos)
  alpha_ecg_heavy <- numeric(tamanho_intervalos)

  sim_ecg_light <- numeric(tamanho_intervalos)
  alpha_ecg_light <- numeric(tamanho_intervalos)

  sim_ecg_dirac <- numeric(tamanho_intervalos)
  alpha_ecg_dirac <- numeric(tamanho_intervalos)

  tipo_ecg_melhor <- character(tamanho_intervalos)

  for (b in seq_len(tamanho_intervalos)) {
    vetor_ecg <- ecg_df$mV[b:(b + janela_tamanho - 1)]

    h.P.sim_heavy <- optim(
      par = c(0.5, 1),
      fn = function(par) Entropy.P(x.var = vetor_ecg, par = par)
    )

    h.P.sim_light <- optim(
      par = c(5, 1),
      fn = function(par) Entropy.P(x.var = vetor_ecg, par = par)
    )

    h.P.sim_dirac <- optim(
      par = c(1000, 1),
      fn = function(par) Entropy.P(x.var = vetor_ecg, par = par)
    )

    sim_ecg_heavy[b] <- h.P.sim_heavy$value
    alpha_ecg_heavy[b] <- h.P.sim_heavy$par[1]

    sim_ecg_light[b] <- h.P.sim_light$value
    alpha_ecg_light[b] <- h.P.sim_light$par[1]

    sim_ecg_dirac[b] <- h.P.sim_dirac$value
    alpha_ecg_dirac[b] <- h.P.sim_dirac$par[1]

    valores <- c(
      heavy = h.P.sim_heavy$value,
      light = h.P.sim_light$value,
      dirac = h.P.sim_dirac$value
    )

    tipo_ecg_melhor[b] <- names(valores)[which.min(valores)]
  }

  blocos_ecg$sim_heavy <- sim_ecg_heavy
  blocos_ecg$alpha_heavy <- alpha_ecg_heavy
  blocos_ecg$sim_light <- sim_ecg_light
  blocos_ecg$alpha_light <- alpha_ecg_light
  blocos_ecg$sim_dirac <- sim_ecg_dirac
  blocos_ecg$alpha_dirac <- alpha_ecg_dirac
  blocos_ecg$tipo_melhor <- tipo_ecg_melhor

  blocos_ecg
}


calcular_entropia_janela_deslizante_ecg_alpha_otimo_com_sobreposicao <- function(
    ecg_df,
    janela_tamanho = 30,
    par_inicial = c(alpha = 5, beta = 1),
    alpha_dirac_limite = 1000
) {
  if (!exists("Entropy.P")) {
    source("EntropyEstimationParetoKernel2.R")
  }

  if (!all(c("time", "mV") %in% names(ecg_df))) {
    stop("O dataframe precisa ter as colunas 'time' e 'mV'.")
  }

  if (janela_tamanho <= 0) {
    stop("O tamanho da janela precisa ser maior que zero.")
  }

  tamanho_intervalos <- nrow(ecg_df) - janela_tamanho + 1

  if (tamanho_intervalos < 1) {
    stop("O dataframe tem menos linhas que o tamanho da janela.")
  }

  blocos_ecg <- data.frame(
    bloco = seq_len(tamanho_intervalos),
    inicio = ecg_df$time[seq_len(tamanho_intervalos)],
    fim = ecg_df$time[seq_len(tamanho_intervalos) + janela_tamanho - 1]
  )

  sim_ecg_otimo <- numeric(tamanho_intervalos)
  alpha_ecg_otimo <- numeric(tamanho_intervalos)
  beta_ecg_otimo <- numeric(tamanho_intervalos)
  tipo_ecg_melhor <- character(tamanho_intervalos)
  convergencia <- integer(tamanho_intervalos)

  classificar_alpha <- function(alpha) {
    if (is.infinite(alpha) || alpha >= alpha_dirac_limite) {
      return("dirac")
    }

    if (alpha <= 2) {
      return("heavy")
    }

    "light"
  }

  for (b in seq_len(tamanho_intervalos)) {
    vetor_ecg <- ecg_df$mV[b:(b + janela_tamanho - 1)]
    penalidade <- 1e100

    objetivo <- function(par) {
      alpha <- par[1]
      beta <- par[2]

      if (!is.finite(alpha) || !is.finite(beta) || alpha <= 0 || beta <= 0) {
        return(penalidade)
      }

      valor <- tryCatch(
        Entropy.P(x.var = vetor_ecg, par = c(alpha, beta)),
        error = function(e) penalidade
      )

      if (!is.finite(valor)) {
        return(penalidade)
      }

      valor
    }

    h.P.sim_otimo <- optim(
      par = par_inicial,
      fn = objetivo,
      method = "L-BFGS-B",
      lower = c(1e-8, 1e-8),
      upper = c(alpha_dirac_limite, Inf)
    )

    sim_ecg_otimo[b] <- h.P.sim_otimo$value
    alpha_ecg_otimo[b] <- h.P.sim_otimo$par[1]
    beta_ecg_otimo[b] <- h.P.sim_otimo$par[2]
    tipo_ecg_melhor[b] <- classificar_alpha(h.P.sim_otimo$par[1])
    convergencia[b] <- h.P.sim_otimo$convergence
  }

  blocos_ecg$sim_otimo <- sim_ecg_otimo
  blocos_ecg$alpha_otimo <- alpha_ecg_otimo
  blocos_ecg$beta_otimo <- beta_ecg_otimo
  blocos_ecg$tipo_melhor <- tipo_ecg_melhor
  blocos_ecg$convergencia <- convergencia

  blocos_ecg
}
