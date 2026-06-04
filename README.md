# Differential Entropy Estimation

Trabalho para a matéria de Tópicos Avançados 2 do Raul, Primeiro Semestre de 2026.

Resumindo:

* Arquivo *main.Rmd*: a simulação principal de monte carlo referente as das tabelas 2, 3, 4 e 5 do artigo original do Raul.

Resultados disso: na pasta *simulaçoes\_14\_05*.

* Arquivo *main\_otimizada.Rmd*: a simulação como antes só que paralelizado com a biblioteca **furr**. Com samples muito grandes, a simulação como no arquivo original nem roda.
* Esses resultados foram então validados no arquivo *analise\_resultados.Rmd*, o que nos indica que os resultados estão dentro do esperado ao comparar com o artigo original.

Resultados disso: na pasta *simulaçoes\_23\_05*.

* Arquivo *bbas3.Rmd*: aplicando entropia nas ações do Banco do Brasil. Isso aqui foi feito muito corrido, provavelmente é necessário revisar.
* Arquivo *main\_otimizada\_pareto.Rmd*: mesma coisa que a main\_otimizada.Rmd só que eu mudei a t-student pela alpha-estável simétrica.

Resultados disso: na pasta *simulacoes\_28\_05*.



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

