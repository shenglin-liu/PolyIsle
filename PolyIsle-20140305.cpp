/*Island model migration,
 k-allele model mutation,
 even-ploidy level,
 hermaphrodite plant*/

#include <iostream>
#include <stdlib.h>		//"srand", "rand"
#include <time.h>		//initialize random seed
#include <random>		//generating random numbers
#include <algorithm>	//"random_shuffle"

int KALLELE=1000;
std::default_random_engine myPRG;
std::uniform_int_distribution<int> myPRDI(1,KALLELE);	//generating random alleles
std::uniform_real_distribution<double> myPRDR(0.0,1.0);	//for migration

/*one locus (an array of genecopies symbolized with integers)*/
class Locus
{
	public:
		Locus();
		Locus(int ploidy, bool randomAssign=true);
		~Locus();
		int equal(Locus & ref);
		Locus & operator=(Locus & ref);
		int & operator[](int temp);

		/*accessor functions*/
		int GetPloidy()const;
		int * GetGenes()const;

		/*model functions*/
		void NewGen(Locus & father, Locus & mother);
			//getting an offsping locus from a father locus and a mother locus
		float GetHo()const;	//calculating Ho

	private:
		int itsPloidy;	//ploidy level
		int * itsGenes;	//array of genecopies
};

/*multiple loci (one individual)*/
class Individual
{
	public:
		Individual();
		Individual(int loci, int ploidy, bool randomAssign=true);
		~Individual();
		int equal(Individual & ref);
		Individual & operator=(Individual & ref);
		Locus & operator[](int temp);

		/*accessor functions*/
		int GetLoci()const;
		int GetPloidy()const;
		Locus * GetGenotype()const;

		/*model functions*/
		void NewGen(Individual & father, Individual & mother);
			//getting an offspring from a father and a mother

	private:
		int itsLoci;	//number of loci
		Locus * itsGenotype;	//array of loci
};

/*multiple individuals (one population)*/
class Population
{
	public:
		Population();
		Population(int popSize, int loci, int ploidy, bool randomAssign=true);
		~Population();
		int equal(Population & ref);
		Population & operator=(Population & ref);
		Individual & operator[](int temp);

		/*accessor functions*/
		int GetPopSize()const;
		int GetLoci()const;
		int GetPloidy()const;
		Individual * GetPop()const;

		/*model functions*/
		float GetHo(int loc)const;	//Ho of specified locus
		float GetHs(int loc)const;	//Hs of specified locus

	private:
		int itsPopSize;		//population size
		Individual * itsPop;	//array of individuals
};

/*multiple populations (one region or meta-population)*/
class Region
{
	public:
		Region();
		Region(int popNum, int popSize, int loci, int ploidy, bool randomAssign=true);
		~Region();
		int equal(Region & ref);
		Region & operator=(Region & ref);
		Population & operator[](int temp);

		/*accessor functions*/
		int GetPopNum()const;
		int GetPopSize()const;
		int GetLoci()const;
		int GetPloidy()const;
		Population * GetRegion()const;

		/*model functions*/
		void Mutate(float mutation);	//mutation
		void NewGen(Region & father, float migration);	//migrate, breed, and drift
		float GetHo(int loc)const;	//Ho of specified locus
		float GetHs(int loc)const;	//Hs of specified locus
		float GetHt(int loc)const;	//Ht of specified locus

	private:
		int itsPopNum;	//number of populations
		Population * itsRegion;	//array of populations
};




/*---------------Implementations for Locus class---------------*/



/*defaut constructor, 
diploid,
initializing genecopies with NULL*/
Locus::Locus()
{
	itsPloidy=2;
	itsGenes=new int[itsPloidy];
	for(int iPloidy=0; iPloidy<itsPloidy; ++iPloidy)
	{
		itsGenes[iPloidy]=0;
	}
}

/*constructing a locus with 
a required ploidy level*/
Locus::Locus(int ploidy, bool randomAssign)
{
	itsPloidy=ploidy;
	itsGenes=new int[itsPloidy];

	if(randomAssign)	//initializing genecopies with random alleles with uniform probability
	{
		for(int iPloidy=0; iPloidy<itsPloidy; ++iPloidy)
		{
			itsGenes[iPloidy]=myPRDI(myPRG);
		}
	}

	else	//initializing genecopies with NULL
	{
		for(int iPloidy=0; iPloidy<itsPloidy; ++iPloidy)
		{
			itsGenes[iPloidy]=0;
		}
	}
}

Locus::~Locus()
{
	delete [] itsGenes;
}

/*deep copy; but the two objects must
be pre-declared*/
int Locus::equal(Locus & ref)
{
	if(this==&ref)
		return 0;
	delete [] itsGenes;
	itsPloidy=ref.GetPloidy();
	itsGenes=new int[itsPloidy];
	for(int iPloidy=0; iPloidy<itsPloidy; ++iPloidy)
	{
		itsGenes[iPloidy]=ref[iPloidy];
	}
	return 1;
}

/*overloading assignment operator;
deep copy without re-allocating the heap;
the two objects must be pre-declared,
have different memory allocations,
and have the same ploidy level*/
Locus & Locus::operator=(Locus & ref)
{
	for(int iPloidy=0; iPloidy<itsPloidy; ++iPloidy)
	{
		itsGenes[iPloidy]=ref[iPloidy];
	}
	return *this;
}

int & Locus::operator[](int temp)
{
	return itsGenes[temp];
}

int Locus::GetPloidy()const
{
	return itsPloidy;
}

int * Locus::GetGenes()const
{
	return itsGenes;
}

/*generating a random gamete from each parent,
combining them to make an offspring locus,
and saving it in an existing object;
the existing object must have a different memory location than the parent objects,
and have the same ploidy level; the parents must have the same ploidy level too, 
but they can be the same object;
NOTE this function changes the orders of the gene arrays of the parents*/
void Locus::NewGen(Locus & father, Locus & mother)
{
	int gametePloidy=itsPloidy/2;
	int iPloidy=0;

	int * sperm=father.GetGenes();
	std::random_shuffle(&sperm[0], &sperm[itsPloidy]);
	for(; iPloidy<gametePloidy; ++iPloidy)
	{
		itsGenes[iPloidy]=sperm[iPloidy];
	}

	int * egg=mother.GetGenes();
	std::random_shuffle(&egg[0], &egg[itsPloidy]);
	for(; iPloidy<itsPloidy; ++iPloidy)
	{
		itsGenes[iPloidy]=egg[iPloidy];
	}
}

/*Ho is calculated through dissimilarity matrix*/
float Locus::GetHo()const
{
	int hetero=0;
	for(int iPloidy=0; iPloidy<itsPloidy; ++iPloidy)
	{
		for(int jPloidy=0; jPloidy<iPloidy; ++jPloidy)
		{
			if(itsGenes[iPloidy]!=itsGenes[jPloidy])
			{
				++hetero;
			}
		}
	}
	return float(hetero)/(itsPloidy*(itsPloidy-1)/2);
}



/*---------------Implementations for Individual class---------------*/



/*default constructor,
diploid, one locus,
initializing genecopies with NULL*/
Individual::Individual()
{
	itsLoci=1;
	itsGenotype=new Locus[itsLoci];
}

/*constructing an individual with
required specifications (number
of loci, ploidy level)*/
Individual::Individual(int loci, int ploidy, bool randomAssign)
{
	itsLoci=loci;
	itsGenotype=new Locus[itsLoci];
	if(randomAssign)	//initializing genecopies with random alleles with uniform probability
	{
		for(int iLoci=0; iLoci<itsLoci; ++iLoci)
		{
			itsGenotype[iLoci].equal(Locus(ploidy));
		}
	}
	else	//initializing genecopies with NULL
	{
		for(int iLoci=0; iLoci<itsLoci; ++iLoci)
		{
			itsGenotype[iLoci].equal(Locus(ploidy, false));
		}
	}
}

Individual::~Individual()
{
	delete [] itsGenotype;
}

/*deep copy; but the two objects must
be pre-declared*/
int Individual::equal(Individual & ref)
{
	if(this==&ref)
		return 0;
	delete [] itsGenotype;
	itsLoci=ref.GetLoci();
	itsGenotype=new Locus[itsLoci];
	for(int iLoci=0; iLoci<itsLoci; ++iLoci)
	{
		itsGenotype[iLoci].equal(ref[iLoci]);
	}
	return 1;
}

/*overloading assignment operator;
deep copy without re-allocating the heap;
the two objects must be pre-declared,
have different memory allocations,
and have the same specifications*/
Individual & Individual::operator=(Individual & ref)
{
	for(int iLoci=0; iLoci<itsLoci; ++iLoci)
	{
		itsGenotype[iLoci]=ref[iLoci];
	}
	return *this;
}

Locus & Individual::operator[](int temp)
{
	return itsGenotype[temp];
}

int Individual::GetLoci()const
{
	return itsLoci;
}

int Individual::GetPloidy()const
{
	return itsGenotype[0].GetPloidy();
}

Locus * Individual::GetGenotype()const
{
	return itsGenotype;
}

/*generating an offspring from parent individuals,
and saving it in an existing object;
the existing object must have a different memory location than the parent objects,
and have the same specifications; the parents must have the same specifications too, 
but they can be the same individual;
NOTE this function changes the orders of
the parents' gene arrays at the locus level*/
void Individual::NewGen(Individual & father, Individual & mother)
{
	for(int iLoci=0; iLoci<itsLoci; ++iLoci)
	{
		itsGenotype[iLoci].NewGen(father[iLoci],mother[iLoci]);
	}
}



/*---------------Implementations for Population class---------------*/



/*default constructor,
diploid, one locus, two individuals
initializing genecopies with NULL*/
Population::Population()
{
	itsPopSize=2;
	itsPop=new Individual[itsPopSize];
}

/*constructing a population with
required specifications (number of individuals,
number of loci, ploidy level)*/
Population::Population(int popSize, int loci, int ploidy, bool randomAssign)
{
	itsPopSize=popSize;
	itsPop=new Individual[itsPopSize];
	if(randomAssign)	//initializing genecopies with random alleles with uniform probability
	{
		for(int iPopSize=0; iPopSize<itsPopSize; ++iPopSize)
		{
			itsPop[iPopSize].equal(Individual(loci, ploidy));
		}
	}
	else	//initializing genecopies with NULL
	{
		for(int iPopSize=0; iPopSize<itsPopSize; ++iPopSize)
		{
			itsPop[iPopSize].equal(Individual(loci, ploidy, false));
		}
	}
}

Population::~Population()
{
	delete [] itsPop;
}

/*deep copy; but the two objects must
be pre-declared*/
int Population::equal(Population & ref)
{
	if(this==&ref)
		return 0;
	delete [] itsPop;
	itsPopSize=ref.GetPopSize();
	itsPop=new Individual[itsPopSize];
	for(int iPopSize=0; iPopSize<itsPopSize; ++iPopSize)
	{
		itsPop[iPopSize].equal(ref[iPopSize]);
	}
	return 1;
}

/*overloading assignment operator;
deep copy without re-allocating the heap;
the two objects must be pre-declared,
have different memory allocations,
and have the same specifications*/
Population & Population::operator=(Population & ref)
{
	for(int iPopSize=0; iPopSize<itsPopSize; ++iPopSize)
	{
		itsPop[iPopSize]=ref[iPopSize];
	}
	return *this;
}

Individual & Population::operator[](int temp)
{
	return itsPop[temp];
}

int Population::GetPopSize()const
{
	return itsPopSize;
}

int Population::GetLoci()const
{
	return itsPop[0].GetLoci();
}

int Population::GetPloidy()const
{
	return itsPop[0].GetPloidy();
}

Individual * Population::GetPop()const
{
	return itsPop;
}

/*at specified locus, calculate average 
Ho of all individuals in the population*/
float Population::GetHo(int loc)const
{
	float ho=0;
	for(int iPopSize=0; iPopSize<itsPopSize; ++iPopSize)
	{
		ho=ho+itsPop[iPopSize][loc].GetHo();
	}
	return ho/itsPopSize;
}

/*Hs is calculated through allele frequency distribution*/
float Population::GetHs(int loc)const
{
	int * freq=new int[KALLELE];
	for(int iFreq=0; iFreq<KALLELE; ++iFreq)
	{
		freq[iFreq]=0;
	}
	int ploidy=this->GetPloidy();
	for(int iPopSize=0; iPopSize<itsPopSize; ++iPopSize)
	{
		for(int iPloidy=0; iPloidy<ploidy; ++iPloidy)
		{
			++freq[itsPop[iPopSize][loc][iPloidy]-1];
		}
	}
	float homo=0; float total=itsPopSize*ploidy; float pi=0;
	for(int iFreq=0; iFreq<KALLELE; ++iFreq)
	{
		pi=freq[iFreq]/total;
		homo=homo+pi*pi;
	}
	delete [] freq;
	return 1-homo;
}



/*---------------Implementations for Region class---------------*/



/*default constructor,
diploid, one locus, two individuals, two populations,
initializing genecopies with NULL*/
Region::Region()
{
	itsPopNum=2;
	itsRegion=new Population[itsPopNum];
}

/*constructing a region with
required specifications (number of populations, 
number of individuals,
number of loci, ploidy level)*/
Region::Region(int popNum, int popSize, int loci, int ploidy, bool randomAssign)
{
	itsPopNum=popNum;
	itsRegion=new Population[itsPopNum];
	if(randomAssign)	//initializing genecopies with random alleles with uniform probability
	{
		for(int iPopNum=0; iPopNum<itsPopNum; ++iPopNum)
		{
			itsRegion[iPopNum].equal(Population(popSize, loci, ploidy));
		}
	}
	else	//initializing genecopies with NULL
	{
		for(int iPopNum=0; iPopNum<itsPopNum; ++iPopNum)
		{
			itsRegion[iPopNum].equal(Population(popSize, loci, ploidy, false));
		}
	}
}

Region::~Region()
{
	delete [] itsRegion;
}

/*deep copy; but the two objects must
be pre-declared*/
int Region::equal(Region & ref)
{
	if(this==&ref)
		return 0;
	delete [] itsRegion;
	itsPopNum=ref.GetPopNum();
	itsRegion=new Population[itsPopNum];
	for(int iPopNum=0; iPopNum<itsPopNum; ++iPopNum)
	{
		itsRegion[iPopNum].equal(ref[iPopNum]);
	}
	return 1;
}

/*overloading assignment operator;
deep copy without re-allocating the heap;
the two objects must be pre-declared,
have different memory allocations,
and have the same specifications*/
Region & Region::operator=(Region & ref)
{
	for(int iPopNum=0; iPopNum<itsPopNum; ++iPopNum)
	{
		itsRegion[iPopNum]=ref[iPopNum];
	}
	return *this;
}

Population & Region::operator[](int temp)
{
	return itsRegion[temp];
}

int Region::GetPopNum()const
{
	return itsPopNum;
}

int Region::GetPopSize()const
{
	return itsRegion[0].GetPopSize();
}

int Region::GetLoci()const
{
	return itsRegion[0].GetLoci();
}

int Region::GetPloidy()const
{
	return itsRegion[0].GetPloidy();
}

Population * Region::GetRegion()const
{
	return itsRegion;
}

/*mutation; decide the number of mutated genecopies
through binomial distribution, randomly choose 
the genecopies across the whole region*/
void Region::Mutate(float mutation)
{
	mutation=mutation*KALLELE/(KALLELE-1);
	int popSize=this->GetPopSize();
	int loci=this->GetLoci();
	int ploidy=this->GetPloidy();
	std::binomial_distribution<int> distribution1(itsPopNum*popSize*loci*ploidy,mutation);
	std::uniform_int_distribution<int> distribution2(0,itsPopNum*popSize*loci*ploidy-1);
	int nMutate=distribution1(myPRG);
	for(int iMutate=0; iMutate<nMutate; ++iMutate)
	{
		int index=distribution2(myPRG);
		int ind1=index%(popSize*loci*ploidy);
		int ind2=ind1%(loci*ploidy);
		int ind3=ind2%ploidy;
		itsRegion[index/(popSize*loci*ploidy)][ind1/(loci*ploidy)][ind2/ploidy][ind3]=myPRDI(myPRG);
	}
}

/*migrate, breed, and drift accoring to the migration rate and the father region,
and save the new region in an existing object;
the existing object must have a different memory location than the father object,
and have the same specifications; 
NOTE this function changes the orders of
the father's gene arrays at the locus level*/
void Region::NewGen(Region & father, float migration)
{
	migration=migration/(1-1/float(itsPopNum));	//correcting for migrating back
	int popSize=father.GetPopSize();
	std::uniform_int_distribution<int> myRandomIndi(0, popSize-1);
	std::uniform_int_distribution<int> myRandomPop(0, itsPopNum-1);

	float pro=0; int jPopNum=0;
	for(int iPopNum=0; iPopNum<itsPopNum; ++iPopNum)
	{
		for(int iPopSize=0; iPopSize<popSize; ++iPopSize)
		{
			pro=myPRDR(myPRG);
			if(pro<migration)
			{
				jPopNum=myRandomPop(myPRG);
				itsRegion[iPopNum][iPopSize].NewGen(father[jPopNum][myRandomIndi(myPRG)], father[jPopNum][myRandomIndi(myPRG)]);
			}
			else
			{
				itsRegion[iPopNum][iPopSize].NewGen(father[iPopNum][myRandomIndi(myPRG)], father[iPopNum][myRandomIndi(myPRG)]);
			}
		}
	}
}

/*at specified locus, calculate average 
Ho of all individuals in all populations*/
float Region::GetHo(int loc)const
{
	float ho=0;
	for(int iPopNum=0; iPopNum<itsPopNum; ++iPopNum)
	{
		ho=ho+itsRegion[iPopNum].GetHo(loc);
	}
	return ho/itsPopNum;
}

/*at specified locus, calculate average 
Hs of all populations*/
float Region::GetHs(int loc)const
{
	float hs=0;
	for(int iPopNum=0; iPopNum<itsPopNum; ++iPopNum)
	{
		hs=hs+itsRegion[iPopNum].GetHs(loc);
	}
	return hs/itsPopNum;
}

/*Ht is calculated through allele frequency distribution*/
float Region::GetHt(int loc)const
{
	int * freq=new int[KALLELE];
	for(int iFreq=0; iFreq<KALLELE; ++iFreq)
	{
		freq[iFreq]=0;
	}
	int popSize=this->GetPopSize();
	int ploidy=this->GetPloidy();
	for(int iPopNum=0; iPopNum<itsPopNum; ++iPopNum)
	{
		for(int iPopSize=0; iPopSize<popSize; ++iPopSize)
		{
			for(int iPloidy=0; iPloidy<ploidy; ++iPloidy)
			{
				++freq[itsRegion[iPopNum][iPopSize][loc][iPloidy]-1];
			}
		}
	}
	float homo=0; float total=itsPopNum*popSize*ploidy; float pi=0;
	for(int iFreq=0; iFreq<KALLELE; ++iFreq)
	{
		pi=freq[iFreq]/total;
		homo=homo+pi*pi;
	}
	delete [] freq;
	return 1-homo;
}



#include <fstream>
#include <string>

using namespace std;

int main()
{
	srand(time(NULL));
	myPRG.seed(rand());	//initializing random number generator in "PopHierarchy.h"

	time_t rawtime;	//timer
	struct tm * timeinfo;

	/*parameters with default values*/
	int popNum=100;
	int popSize=1000;
	int loci=10;
	float mutation=0.00001;
	
	/*lengths to be defined*/
	int burnin=-1;
	int generation=-1;
	int echo=-1;	//how often to echo to screen
	int save=-1;	//how often to calculate and save the heterozygosities

	/*variables to be defined*/
	float migration=-1;	//migration rate during real simulation
	float migrationB=1;	//migration rate during burn-in, default as panmictic
	int ploidy=-1;

	/*initializing output files*/
	ofstream fProfile("profile.txt");	//experimental parameters and time of start/burnin-complete/simulation-complete
	ofstream fBurnin("burnin.txt");	//output of burnin
	ofstream fSim("simulation.txt");	//output of real simulation
	ofstream fState;	//output a region object

	/*initializing the experiment (using commandlines)*/
	int cBurnin=0, cGen=0, cEcho=0, cSave=0, cMig=0, cPloidy=0;	//checkpoints
	int check;	//check whether the experiment is properly initialized
	int append=0;	//continue an existing experiment(1), or not(0)
	string command;
	cout<<"Please set up your experiment using commandlines.\n";
	for(;;)
	{
		cout<<">>";
		cin>>command;
		if(command=="popNum")
		{
			cin>>popNum;
			continue;
		}
		if(command=="popSize")
		{
			cin>>popSize;
			continue;
		}
		if(command=="loci")
		{
			cin>>loci;
			continue;
		}
		if(command=="mutation")
		{
			cin>>mutation;
			continue;
		}
		if(command=="burnin")
		{
			cin>>burnin;
			cBurnin=1;
			continue;
		}
		if(command=="generation")
		{
			cin>>generation;
			cGen=1;
			continue;
		}
		if(command=="echo")
		{
			cin>>echo;
			cEcho=1;
			continue;
		}
		if(command=="save")
		{
			cin>>save;
			cSave=1;
			continue;
		}
		if(command=="migration")
		{
			cin>>migration;
			cMig=1;
			continue;
		}
		if(command=="migrationB")
		{
			cin>>migrationB;
			continue;
		}
		if(command=="ploidy")
		{
			cin>>ploidy;
			cPloidy=1;
			continue;
		}
		if(command=="append")
		{
			cin>>append;
			continue;
		}
		if(command=="summary")
		{
			cout<<"popNum=\t\t"<<popNum<<endl;
			cout<<"popSize=\t"<<popSize<<endl;
			cout<<"loci=\t\t"<<loci<<endl;
			cout<<"mutation=\t"<<mutation<<endl;
			cout<<"burnin=\t\t"<<burnin<<endl;
			cout<<"generation=\t"<<generation<<endl;
			cout<<"echo=\t\t"<<echo<<endl;
			cout<<"save=\t\t"<<save<<endl;
			cout<<"migration=\t"<<migration<<endl;
			cout<<"migrationB=\t"<<migrationB<<endl;
			cout<<"ploidy=\t\t"<<ploidy<<endl;
			cout<<"append=\t\t"<<append<<endl;
			continue;
		}
		if(command=="start")
		{
			check=cBurnin*cGen*cEcho*cSave*cMig*cPloidy;
			if(check==0)
			{
				cout<<"Experimental setting incomplete! Please check using 'summary' command!\n";
				continue;
			}
			else
			{
				string ready;
				cout<<"Are you ready to start?[Y/N]\n";
				cin>>ready;
				if(ready=="Y")
				{
					cin.ignore(50,'\n');
					break;
				}
				if(ready=="N")
				{
					cin.ignore(50,'\n');
					continue;
				}
				else
				{
					cout<<"Wrong answer!\n";
					cin.ignore(50,'\n');
				}
			}
		}
		else
		{
			cout<<"Wrong command!\n";
			cin.ignore(50,'\n');
		}
	}

	/*recording the settings and the starting time of the experiment*/
	fProfile<<"popNum=\t\t"<<popNum<<endl;
	fProfile<<"popSize=\t"<<popSize<<endl;
	fProfile<<"loci=\t\t"<<loci<<endl;
	fProfile<<"mutation=\t"<<mutation<<endl;
	fProfile<<"burnin=\t\t"<<burnin<<endl;
	fProfile<<"generation=\t"<<generation<<endl;
	fProfile<<"echo=\t\t"<<echo<<endl;
	fProfile<<"save=\t\t"<<save<<endl;
	fProfile<<"migration=\t"<<migration<<endl;
	fProfile<<"migrationB=\t"<<migrationB<<endl;
	fProfile<<"ploidy=\t\t"<<ploidy<<endl;
	fProfile<<"append=\t\t"<<append<<endl;
	time(&rawtime);
	timeinfo=localtime(&rawtime);
	fProfile<<"\nExperiment started at: "<<asctime(timeinfo)<<endl;

	/*start a region object(population structure)*/
	Region region(popNum,popSize,loci,ploidy);		//real object to collect data from
	Region buffer(popNum,popSize,loci,ploidy,false);	//buffer object for migration and reproducing
	if(append!=0)	//importing the region, if you have chosen to continue an existing experiment
	{
		ifstream stateApp("state.txt");
		for(int i=0; i<popNum; ++i)
		{
			for(int j=0; j<popSize; ++j)
			{
				for(int k=0; k<loci; ++k)
				{
					for(int l=0; l<ploidy; ++l)
					{
						stateApp>>region[i][j][k][l];
					}
				}
			}
		}
		stateApp.close();
		cout<<"Region object import complete!"<<endl;
	}
	
	/*heterozygosities*/
	float ho=0, hs=0, ht=0;

	/*burnin*/
	int actBurnin=burnin/2;
	for(int i=1; i<=actBurnin; ++i)
	{
		region.Mutate(mutation);	//mutate
		buffer.NewGen(region, migrationB);	//breed and drift

		if(((i*2-1)%echo)==0)	//echo to the screen
		{
			cout<<"burnin "<<i*2-1<<endl;
		}

		if(((i*2-1)%save)==0)	//calculate heterozygosities and save
		{
			ho=0; hs=0; ht=0;
			for(int j=0; j<loci; ++j)
			{
				ho=ho+buffer.GetHo(j);
				hs=hs+buffer.GetHs(j);
				ht=ht+buffer.GetHt(j);
			}
			ho=ho/loci;
			hs=hs/loci;
			ht=ht/loci;
			fBurnin<<ho<<"\t"<<hs<<"\t"<<ht<<"\t"<<endl;
		}
		
		buffer.Mutate(mutation);
		region.NewGen(buffer, migrationB);

		if(((i*2)%echo)==0)
		{
			cout<<"burnin "<<i*2<<endl;
		}

		if(((i*2)%save)==0)
		{
			ho=0; hs=0; ht=0;
			for(int j=0; j<loci; ++j)
			{
				ho=ho+region.GetHo(j);
				hs=hs+region.GetHs(j);
				ht=ht+region.GetHt(j);
			}
			ho=ho/loci;
			hs=hs/loci;
			ht=ht/loci;
			fBurnin<<ho<<"\t"<<hs<<"\t"<<ht<<"\t"<<endl;
		}
	}
	fBurnin.close();	//burnin complete!
	time(&rawtime);	//record time
	timeinfo=localtime(&rawtime);
	fProfile<<"\nBurnin completed at: "<<asctime(timeinfo)<<endl;

	/*real simulation*/
	int actSimu=generation/2;
	for(int i=1; i<=actSimu; ++i)
	{
		region.Mutate(mutation);	//mutate
		buffer.NewGen(region, migration);	//migrate, breed, and drift
		
		if(((i*2-1)%echo)==0)	//echo to the screen
		{
			cout<<"simulation "<<i*2-1<<endl;
		}

		if(((i*2-1)%save)==0)	//calculate heterozygosities and save
		{
			ho=0; hs=0; ht=0;
			for(int j=0; j<loci; ++j)
			{
				ho=ho+buffer.GetHo(j);
				hs=hs+buffer.GetHs(j);
				ht=ht+buffer.GetHt(j);
			}
			ho=ho/loci;
			hs=hs/loci;
			ht=ht/loci;
			fSim<<ho<<"\t"<<hs<<"\t"<<ht<<"\t"<<endl;
		}
		
		buffer.Mutate(mutation);
		region.NewGen(buffer, migration);

		if(((i*2)%echo)==0)
		{
			cout<<"simulation "<<i*2<<endl;
		}

		if(((i*2)%save)==0)
		{
			ho=0; hs=0; ht=0;
			for(int j=0; j<loci; ++j)
			{
				ho=ho+region.GetHo(j);
				hs=hs+region.GetHs(j);
				ht=ht+region.GetHt(j);
			}
			ho=ho/loci;
			hs=hs/loci;
			ht=ht/loci;
			fSim<<ho<<"\t"<<hs<<"\t"<<ht<<"\t"<<endl;
		}
	}
	fSim.close();	//simulation complete!
	time(&rawtime);	//record time
	timeinfo=localtime(&rawtime);
	fProfile<<"\nSimulation completed at: "<<asctime(timeinfo)<<endl;
	fProfile.close();

	fState.open("state.txt", fstream::trunc);	//save the final state of the region object
	for(int j=0; j<popNum; ++j)
	{
		for(int k=0; k<popSize; ++k)
		{
			for(int l=0; l<loci; ++l)
			{
				for(int m=0; m<ploidy; ++m)
				{
					fState<<region[j][k][l][m]<<"\t";
				}
				fState<<endl;
			}
		}
	}
	fState.close();

	cout<<"Experiment is complete!\nInput any character and press enter to terminate the application."<<endl;

	char response;
	cin>>response;

	return 0;
}