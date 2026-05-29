# Differential Entropy Estimation

Trabalho para a matéria de Tópicos Avançados 2 do Raul, Primeiro Semestre de 2026.

Resumindo:

-   Arquivo *main.Rmd*: a simulação principal de monte carlo referente as das tabelas 2, 3, 4 e 5 do artigo original do Raul.

Resultados disso: na pasta *simulaçoes_14_05*.

-   Arquivo *main_otimizada.Rmd*: a simulação como antes só que paralelizado com a biblioteca **furr**. Com samples muito grandes, a simulação como no arquivo original nem roda.

Resultados disso: na pasta *simulaçoes_23_05*.

-   Arquivo *bbas3.Rmd*: aplicando entropia nas ações do Banco do Brasil. Isso aqui foi feito muito corrido, provavelmente é necessário revisar.
 
-   Arquivo *main_otimizada_pareto.Rmd*: mesma coisa que a main_otimizada.Rmd só que eu mudei a t-student pela alpha-estável simétrica.

Resultados disso: na pasta *simulacoes_28_05*.