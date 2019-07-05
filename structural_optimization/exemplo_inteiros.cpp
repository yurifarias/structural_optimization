/* ----------------------------------------------------------------------------
  exemplo_inteiros.cpp
  yurifarias 04/07/19

 DESCRIÇÃO:
	Programa que busca os indices do pilares que recebem as menores cargas

---------------------------------------------------------------------------- */
#include <stdio.h>
#include <ga/ga.h>
#include <ga/std_stream.h>
#include <algorithm>
#include <vector>

#define cout STD_COUT

using namespace std;

float objective(GAGenome&);
int carga_pilar[] = { 11 * 50,	// 0 
					  20 * 50,	// 1
					  7 * 50,	// 2 ---
					  3 * 50,	// 3 ---
					  19 * 50,	// 4
					  8 * 50,	// 5 ---
					  25 * 50,	// 6
					  10 * 50,	// 7
					  1 * 50,	// 8 ---
					  15 * 50,	// 9
					  21 * 50,	// 10
					  14 * 50,	// 11
					  9 * 50,	// 12
					  23 * 50,	// 13
					  4 * 50,	// 14 ---
					  12 * 50,	// 15
					  2 * 50,	// 16 ---
					  18 * 50,	// 17
					  5 * 50,	// 18 ---
					  24 * 50,	// 19
					  16 * 50,	// 20
					  22 * 50,	// 21
					  13 * 50,	// 22
					  6 * 50,	// 23 ---
					  17 * 50	// 24
					}; // valores de 1 a 25 embaralhados

int main(int argc, char** argv)
{
	cout << "Programa de Pilares\n" << endl;
	cout << "Dado o vetor de cargas de pilares: \n\n" << endl;
	for (int i = 0; i < 25; i++)
	{
		cout << i << ": " << carga_pilar[i] << endl;
	}
	cout << "\n\nO programa tenta achar os oito pilares que" << endl;
	cout << "recebem as menores cargas por seus indices" << endl;
	cout << "\n\n"; cout.flush();

	// See if we've been given a seed to use (for testing purposes).  When you
	// specify a random seed, the evolution will be exactly the same each time
	// you use that seed number.

	unsigned int seed = 0;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i++], "seed") == 0) {
			seed = atoi(argv[i]);
		}
	}

	// Declare variables for the GA parameters and set them to some default values.

	int length = 8;
	int popsize = 100;
	int ngen = 400;
	float pmut = 0.20;
	float pcross = 0.20;
	int mm = -1;

	// crio um array de tamanho 25 (tamanho do array de cargas acima)
	int aset[25];

	// preencho o array com os indices do array de cargas
	for (int i = 0; i <= 24; i++) {
		aset[i] = i;
	}

	GAAlleleSet<int> allele(25, aset);
	GA1DArrayAlleleGenome<int> genome(length, allele, objective);

	GASimpleGA ga(genome);
	GASigmaTruncationScaling scaling;
	ga.populationSize(popsize);
	ga.nGenerations(ngen);
	ga.pMutation(pmut);
	ga.pCrossover(pcross);
	ga.scaling(scaling);
	ga.minimaxi(mm);
	ga.scoreFilename("bog.dat");
	ga.scoreFrequency(10);
	ga.flushFrequency(50);
	ga.evolve(seed);

	// Dump the results of the GA to the screen.
	cout << "O programa encontrou os indices:\n" << endl;
	cout << ga.statistics().bestIndividual();
	cout << "\n\n"; cout.flush();

	return 0;
}

float objective(GAGenome& c)
{
	GA1DArrayAlleleGenome<int>& genome = (GA1DArrayAlleleGenome<int>&)c;


	// score (vamos achar o valor minimo)
	int cargas = { 0 };
	int multiplicador = { 1 };

	// inicializo um vetor para ser preenchido sem repetição de valores
	vector <int> indices_pilares;

	// adiciono ao vetor o valor do primeiro gene
	indices_pilares.push_back(genome.gene(0));

	// somo as cargas considerando o multiplicador
	cargas += carga_pilar[genome.gene(0)] * multiplicador;

	// loop para adicionar mais 7 indices de pilares (sem repetir valor)
	for (int i = 1; i < 8; i++)
	{
		// gero um valor aleatorio entre 0 e 24
		int indice_para_adicionar = genome.gene(i);

		// checo se ele ja pertence ao vetor de indices de pilares
		auto valor_diferente = find(begin(indices_pilares), end(indices_pilares), indice_para_adicionar);

		if (valor_diferente != end(indices_pilares)) {
			// se o valor ja existir
			multiplicador = 100;
		}
		else {
			// se nao existir
			multiplicador = 1;
		}

		// adiciono o valor pro final do vetor
		indices_pilares.push_back(indice_para_adicionar);

		// somo as cargas considerando o multiplicador
		cargas += carga_pilar[genome.gene(i)] * multiplicador;
	}

	return cargas;
}
