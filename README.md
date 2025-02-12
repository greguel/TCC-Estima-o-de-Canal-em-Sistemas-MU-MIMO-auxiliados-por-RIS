# TCC-Estima-o-de-Canal-em-Sistemas-MU-MIMO-auxiliados-por-RIS
Código fonte usado no meu TCC comparando métodos de estimação de canal em sistemas MU-MIMO auxiliados por RIS.

O código é composto de quatro partes principais de algoritmo e alguns objetos auxiliares. Resumidamente, os parâmetros do
canal são definidos em um objeto do tipo Channel, e os objetos SimulationDFT e Simulation3P são responsáveis por fazer uma 
única vez os respecitvos métodos de estimação. O objeto do tipo MonteCarlo é responsável por rodar a simulação o número de 
vezes desejado para cada parâmetro de entrada, mudar o parâmetro, e depois plotar os resultados. O arquivo New_RIS é o onde
os parâmetros são definidos e chamados. Alguns dos parâmetros possuem valor padrão na declaração de Channel.
