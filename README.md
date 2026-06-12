# Differential Entropy Estimation

Trabalho para a matéria de Tópicos Avançados 2 do Raul, Primeiro Semestre de 2026.

# Arquivos:

## Exploração

* Arquivo *main.Rmd*: a simulação principal de monte carlo referente as das tabelas 2, 3, 4 e 5 do artigo original do Raul.

Resultados disso: na pasta *simulaçoes\_14\_05*.

* Arquivo *main\_otimizada.Rmd*: a simulação como antes só que paralelizado com a biblioteca **furr**. Com samples muito grandes, a simulação como no arquivo original nem roda.
* Esses resultados foram então validados no arquivo *analise\_resultados.Rmd*, o que nos indica que os resultados estão dentro do esperado ao comparar com o artigo original.

Resultados disso: na pasta *simulaçoes\_23\_05*.

* Arquivo *bbas3.Rmd*: aplicando entropia nas ações do Banco do Brasil. Com isso dá pra tentar reproduzir as Figuras 3 e 4 do Artigo original.

Resultados: pasta *simulacoes_04_06*.

* Arquivo *main\_otimizada\_pareto.Rmd*: mesma coisa que a main\_otimizada.Rmd só que eu mudei a t-student pela alpha-estável simétrica.

Resultados disso: na pasta *simulacoes\_28\_05*.

Até aqui são os resultados das simulações de Monte Carlo da parte teórica do trabalho.

## O que acabou indo pro poster, entrega da matéria em junho de 2026

* Arquivo *bbas3.Rmd*: a análise final que foi pro poster referente a análise de bbas3 com janela fixa e gráficos.

* Arquivo *entropia_paciente_saudavel.Rmd* e *entropia_paciente_saudavel.Rmd*: a contagem de tipos de caudas (leve e pesada) são as que foram pro poster, referente a paciente saudável e doente (no caso um paciente considerado muito doente e outro baseline simples).

* Arquivo *entropia_paciente_cesar.Rmd* usa a base ecgs_por_paciente_cesar com um exemplo de paciente com infarto, com isquemia e um saudável para gerar análises de variância pelp exponencial da entropia com o objetivo de tentar agrupar pacientes.

* Dentro da pasta *dados_aplicacao* tem 2 scripts: um que gera um .rds da base inteira pronta para consumo e um que le os diagnósticos.

* O arquivo *entropia_por_paciente* gera a entropia de um paciente e salva num arquivo chamado resultado_entropia_paciente.rds sempre dando append.

Usei esse script para gerar a base inteira que pode ser baixada como um .rds aqui nesse [link](https://unbbr-my.sharepoint.com/:u:/g/personal/261102693_aluno_unb_br/IQBH6h6WJq_HTJL4ZxYzFlO3ARfLx_D8DbvuNlWIgTaCA1U?e=2axQor)

## Dados utilizados para aplicação

https://physionet.org/content/ptb-xl/1.0.3/ -- download disponível no final da página

* Arquivo *dados_aplicacao/read_ecg_data.R*: lê dados para um id de paciente, com a série de todas as derivações (perspectivas sobre o miocárdio) e possibilidade de montar o gráfico do ECG
* Arquivo *dados_aplicacao/read_patient_diagnostics.R*: gera uma tabela com dados de interesse para filtrarmos casos. Queremos confidence = 100 e comparar scp_ecg_statement_desc == "normal ECG" com "ischemic ST-T changes" ou "(complete) left bundle branch block"


## Citações necessárias para a fonte de dados

@article{PhysioNet-ptb-xl-1.0.3,
  author = {Wagner, Patrick and Strodthoff, Nils and Bousseljot, Ralf-Dieter and Samek, Wojciech and Schaeffter, Tobias},
  title = {{PTB-XL, a large publicly available electrocardiography dataset}},
  journal = {{PhysioNet}},
  year = {2022},
  month = nov,
  note = {Version 1.0.3},
  doi = {10.13026/kfzx-aw45},
  url = {https://doi.org/10.13026/kfzx-aw45}
}

Additionally, please cite the original publication:

Wagner, P., Strodthoff, N., Bousseljot, R.-D., Kreiseler, D., Lunze, F.I., Samek, W., Schaeffter, T. (2020), PTB-XL: A Large Publicly Available ECG Dataset. Scientific Data. https://doi.org/10.1038/s41597-020-0495-6

Please include the standard citation for PhysioNet: (show more options)

Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... \& Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation \[Online]. 101 (23), pp. e215–e220. RRID:SCR\_007345. 

