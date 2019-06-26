/* ----------------------------------------------------------------------------
  exemplo1.cpp
  yurifarias 22/06/19
  atualização1 25/06/19
  
 DESCRIÇÃO:
	Programa que realiza a otimização de uma viga biapoiada de concreto armado,
	cujos dados passados são:

		o momento atuante Md  (kN * cm);
		a tensão resistente do concreto fcd (kN / cm^2);
		a tensão resistente do aço sigmaSd (kN / cm^2);
		o coeficiente betaX (x / d);
		a altura útil d (cm);
		o cobrimento (cm);
		o coeficiente para o aço alphaA (R$/cm^3);
		o coeficiente para o concreto alphaC (R$/cm^3);
		o coeficiente para o forma alphaF (R$/cm^2)

---------------------------------------------------------------------------- */
#include <stdio.h>
#include <ga/ga.h>
#include <ga/std_stream.h>

#define cout STD_COUT

float objective(GAGenome&);
float calcularMomento(float, float);
float calcularAreaDeAco(float, float);

float momAtuante = 10000; // kN * cm
float sigmaSd = 50000 / 1.15; // 50000 kN / cm^2
float cobrimento = 2; // cm
float betaX = 0.45;

int
main(int argc, char** argv)
{
	cout << "Viga biapoiada\n";
	cout << "O programa tenta achar o melhor valor para a altura e a base da secao\n";
	cout << "de forma que minimize a funcao custo:\n";
	cout << "fCusto = alphaAco * areaDeAco + \n";
	cout << "         alphaConcreto * (base * altura - areaDeAco)\n";
	cout << "         alphaForma * 2 * (base + altura)\n";
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

	int popsize = 30;
	int ngen = 100;
	float pmut = 0.01;
	float pcross = 0.6;
	int mm = -1;

	// Create a phenotype for two variables.  The number of bits you can use to
	// represent any number is limited by the type of computer you are using.  In
	// this case, we use 16 bits to represent a floating point number whose value
	// can range from -5 to 5, inclusive.  The bounds on x1 and x2 can be applied
	// here and/or in the objective function.

	// Aqui crio o fenótipo e adiciono minhas variáveis,
	// sendo o primeiro a base da viga e o segundo a altura.

	GABin2DecPhenotype map;
	map.add(8, 10, 25); // Limitei a base da viga para estar entre 10cm e 25cm
	map.add(8, 30, 100); // Limitei a altura da viga para estar entre 30cm e 100cm

	// Create the template genome using the phenotype map we just made.

	// Aqui é o conversor de bits para decimais

	GABin2DecGenome genome(map, objective);

	// Now create the GA using the genome and run it.  We'll use sigma truncation
	// scaling so that we can handle negative objective scores.

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

	genome = ga.statistics().bestIndividual();
	float bestMom = calcularMomento(genome.phenotype(0), genome.phenotype(1));
	float bestAs = calcularAreaDeAco(bestMom, genome.phenotype(1));
	cout << "Os melhores valores para a secao: (";
	cout << genome.phenotype(0) << ", " << genome.phenotype(1) + cobrimento << ")\n";
	cout << "area de aco encontrada: " << bestAs << " cm^2.\n\n";
	cout << "melhor solucao salva em: '" << ga.scoreFilename() << "'\n";

	return 0;
}


// Minha função vai ser uma função de custo em função da área da seção
//
//				custo = alphaAço * areaDeAço + 
//						alphaConcreto * (base * altura - areaDeAço) + 
//						alphaForma * 2 * (base + altura)
//
/*
			Md = 0.68 * bw * betaX * d^2 * fcd * (1 - 0.4 * betaX)
onde:

	momAtuante = momento de cálculo atuante (kN * cm)
	momReativo = momento de cálculo resistente (kN * cm)
	fcd = resitência do concreto (kN / cm^2)
	sigmaSd = resitência do aço (kN / cm^2)
	bw = base da seção (cm)
	d = altura útil (cm)
	betaX = x / d
	x = altura da linha neutra (cm)
	h = altura total (cm)
	h = d + c
	c = cobrimento (cm)

			As = Md / (sigmaS * (1 - 0.4 * betaX))
*/

float
objective(GAGenome& c)
{
	GABin2DecGenome& genome = (GABin2DecGenome&)c;
	float bw = genome.phenotype(0);
	float d = genome.phenotype(1);

	float momReativo = calcularMomento(bw, d);

	float as = calcularAreaDeAco(momReativo, d);

	// aqui testo se o momento reativo é menor que o momento atuante
	// caso seja, estouro os valores de bw e h para que resultem um custo alto
	// e não sejam selecionadas como possível solução.
	if (momReativo < momAtuante)
	{
		bw = bw * 100;
		d = d * 100;
	}

	float custo; // R$/cm
	float alphaC = 200; // sei lá, 100 R$/cm^3
	float alphaA = 5000; // R$/cm^3
	float alphaF = 50; // R$/cm^2

	// funcao custo
	custo = alphaA * as + alphaC * (bw * (d + cobrimento) - as) + alphaF * 2 * (bw + d);

	return custo;
}

float
calcularMomento(float bw, float d) 
{
	float momReativo;

	float fcd = 2.0 / 1.4; // kN / cm^2
	float h = d + cobrimento; // cm

	// função tirada da apostila do paulo bastos (Eq. 20)
	momReativo = 0.68 * bw * betaX * d * d * fcd * (1 - 0.4 * betaX);

	return momReativo;
}

float
calcularAreaDeAco(float momReativo, float d)
{
	float as;

	// área de aço da armadura tracionada
	as = momReativo / (sigmaSd * (1 - 0.4 * betaX) * d);

	return as;
}
