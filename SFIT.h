#pragma once
#ifndef SFIT_H
#define SFIT_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <ctime>

using namespace std;

struct IDnode
{
	short ID;
	IDnode * next;
};

enum LayerType
{
	FIRST, INTER, LAST
};
struct response
{
	float fuzzyValue;
	char classID;
	float probability;
	string s;
};

struct layer
{
	unsigned char ID;
	unsigned int sizeX, sizeY;
	int **index;
	unsigned char **fuzzy;
	unsigned short paramNo0;
	unsigned short paramNo1;
	unsigned int contents; //number of registered combinations in a layer
							 //LayerType lType;
	double a; //linear transforming coefficients
	int b;
	layer * next;

};

struct numberNode
{
	int data;
	numberNode * next;
};

class SFIT
{
private:
	layer *first;

	layer *last;

	int NoC;
	int counter;
	double * currentSample;
	string * patternList;
	float * A;
	int * B;
	int * D;
	int NoE; //number of class entries
	IDnode **items;
	numberNode *numbersFirst;
	numberNode *numbersLast;
	int numberCounter;
	int * intervals;
public:

	SFIT()
	{
		first = NULL; counter = 0; last = NULL; A = NULL; B = NULL; currentSample = NULL; patternList = NULL;//patternList = new PatternList();
		items = NULL; numbersFirst = NULL; numbersLast = NULL; numberCounter = 0; intervals = NULL;
	}
	void initializeLayer(layer * l, int num0, int num1)
	{
		l->sizeX = num0;
		l->sizeY = num1;

		l->index = new int*[num0];
		l->fuzzy = new unsigned char*[num0];
		for (int i = 0; i<num0; i++)
		{
			l->index[i] = new int[num1];
			l->fuzzy[i] = new unsigned char[num1];

			for (int j = 0; j<num1; j++)
			{
				l->index[i][j] = 0;
				l->fuzzy[i][j] = 0;
			}
		}
	}

	void initializeLayer(layer * l, int num0, int num1, int init)
	{
		l->sizeX = num0;
		l->sizeY = num1;

		l->index = new int*[num0];
		l->fuzzy = new unsigned char*[num0];
		for (int i = 0; i<num0; i++)
		{
			l->index[i] = new int[num1];
			l->fuzzy[i] = new unsigned char[num1];

			for (int j = 0; j<num1; j++)
			{
				l->index[i][j] = init;
				l->fuzzy[i][j] = 0;
			}
		}

	}
	void appendFirst(int num0, short param0, int init)
	{
		if (first != NULL) return;

		first = new layer;
		first->ID = 0;
		initializeLayer(first, num0, 1, -1);

		first->paramNo0 = param0;
		first->contents = 0;
		first->next = NULL;
		last = first;
		counter += 1;
	}
	void appendFirst(int num0, short param0)
	{
		cout << "\n gotten paramNo:" << param0;
		if (first != NULL) return;

		first = new layer;
		first->ID = 0;
		initializeLayer(first, num0, 1);
		//first->classes = NULL;
		//first->classCounter = NULL;
		//first->mostClass = NULL;
		//first->lType = FIRST;

		first->paramNo0 = param0;
		first->contents = 0;
		first->next = NULL;
		last = first;
		counter += 1;
	}
	void appendInter(int num0, int num1, short param)
	{
		layer *t;
		t = new layer;

		initializeLayer(t, num0, num1);
		t->paramNo0 = param;
		t->paramNo1 = -1;
		t->contents = 0;
		t->next = NULL;
		last->next = t;
		last = t;
		t->ID = counter;
		counter += 1;
	}
	void appendInter(int num0, int num1, short param, int init)
	{
		layer *t;
		t = new layer;

		initializeLayer(t, num0, num1, -1);
		t->paramNo0 = param;
		t->paramNo1 = -1;
		t->contents = 0;
		t->next = NULL;
		last->next = t;
		last = t;
		t->ID = counter;
		counter += 1;
	}

	void runDiagnostics()
	{
		layer *q = first;
		cout << "\n DIAGNOSTICS";
		cout << "\n-------------";
		while (q)
		{
			cout << "\n > Layer#" << (int)q->ID << ", size=" << q->sizeX << " x " << q->sizeY << "    -    contents:" << q->contents;//<<"    -    contents:"<<q->;
			q = q->next;
		}
	}
	void listContent()
	{
		layer *q = first;
		cout << "\n Contents";
		cout << "\n-------------";
		while (q)
		{
			cout << "\n Layer #" << q->ID;
			for (int i = 0; i<q->sizeX; i++)
			{
				cout << "\n";
				for (int j = 0; j<q->sizeY; j++)
				{
					cout << " " << q->index[i][j];
				}
			}
			for (int i = 0; i<q->sizeX; i++)
			{
				cout << "\n";
				for (int j = 0; j<q->sizeY; j++)
				{
					cout << " " << (int)q->fuzzy[i][j];
				}
			}
			q = q->next;
		}
	}	
	void listContentPartial()
	{
		layer *q = first;
		cout << "\n Contents";
		cout << "\n-------------";
		while (q)
		{
			cout << "\n Layer #" << q->ID;
			for (int i = 0; i<10; i++)
			{
				cout << "\n";
				for (int j = 0; j<q->sizeY; j++)
				{
					cout << " " << q->index[i][j];
				}
			}
			for (int i = q->sizeX-10; i<q->sizeX; i++)
			{
				cout << "\n";
				for (int j = 0; j<q->sizeY; j++)
				{
					cout << " " << q->index[i][j];
				}
			}
			for (int i = 0; i<10; i++)
			{
				cout << "\n";
				for (int j = 0; j<q->sizeY; j++)
				{
					cout << " " << (int)q->fuzzy[i][j];
				}
			}
			for (int i = q->sizeX-10; i<q->sizeX; i++)
			{
				cout << "\n";
				for (int j = 0; j<q->sizeY; j++)
				{
					cout << " " << (int)q->fuzzy[i][j];
				}
			}
			q = q->next;
		}
	}


	void logContent(char * filename)
	{
		ofstream myfile;
		myfile.open(filename);
		layer *q = first;
		myfile << "\n Contents";
		myfile << "\n-------------";
		myfile << "\n Layer #" << (int)q->ID;
		myfile << "\n  - index:";
		for (int j = 0; j<q->sizeY; j++)
		{
			myfile << " " << q->index[j][0];
		}

		myfile << "\n  - Fuzzy:";
		for (int j = 0; j<q->sizeY; j++)
		{
			myfile << " " << (int)q->fuzzy[j][0];
		}
		q = q->next;
		while (q)
		{
			myfile << "\n Layer #" << (int)q->ID;
			myfile << "\n  - index:";
			for (int i = 0; i<q->contents; i++)
			{
				myfile << "\n";
				for (int j = 0; j<q->sizeY; j++)
				{
					myfile << " " << q->index[i][j];
				}
			}
			myfile << "\n  - Fuzzy:";
			for (int i = 0; i<q->contents; i++)
			{
				myfile << "\n";
				for (int j = 0; j<q->sizeY; j++)
				{
					myfile << " " << (int)q->fuzzy[i][j];
				}
			}
			q = q->next;
		}
		myfile.close();
	}
	void calculateAndSetFuzzy1D(layer * L, int row, int center, int RO) // calculate the RO-sized area around (row, center) of L
	{
		for (int i = max(0, center - RO); i<min(center + RO + 1, (int)L->sizeY); i++)
		{
			float aa = (abs(center - i));
			aa /= RO;
			aa = 1 - aa;
			aa *= 100;
			if (L->fuzzy[row][i] < aa)
			{
				L->fuzzy[row][i] = aa;
				L->index[row][i] = L->index[row][center];
			}
		}
		L->fuzzy[row][center] = 100;
	}
	void calculateAndSetFuzzy2D(layer * L, int row, int center, short RO0, short RO1) // calculate the RO-sized area around (row, center) of L
	{
		float a;
		for (int i = max(0, row - RO0); i<min(row + RO0 + 1, (int)L->sizeX); i++)
			for (int j = max(0, center - RO1); j<min(center + RO1 + 1, (int)L->sizeY); j++)
			{
				a = 1 - sqrt(static_cast<float>((row - i)*(row - i) + (center - j)*(center - j))) / ((RO0 + RO1) / 2);
				a *= 100;
				if (L->fuzzy[i][j] < a)
				{
					L->fuzzy[i][j] = a;
					L->index[i][j] = L->index[row][center];
				}
			}
		L->fuzzy[row][center] == 100;
	}

	void train(char* trainingFile, int display, float RP)
	//only needs the scales, automatically calculates B and D, fuzzy width percentage RP
	{
		fstream myfile;
		myfile.open(trainingFile, ios::in);

		int P, N, RO;
		RO = 0;

		myfile >> P;
		myfile >> N;

		int aa;
		double ff;

		A = new float[N];

		for (int i = 0; i < N; i++)
		{
			myfile >> ff;
			A[i] = ff;
		}
		if (display)
		{
			cout << "\n ";
			for (int i = 0; i < N; i++)
			{
				cout << " " << A[i];
			}
		}

		if (display) cout << "\n ";

		int ** trainingValues = new int*[P];
		for (int i = 0; i < P; i++)
		{
			trainingValues[i] = new int[N];
		}
		for (int p = 0; p < P; p++)
		{
			//look at the first layer
			for (int i = 0; i < N; i++)
			{
				myfile >> ff;
				trainingValues[p][i] = round(ff*A[i]);
			}
		}
		myfile.close();

		clock_t start, end;
		double elapsed_sec;
		start = clock();

		if (display)
		{
			for (int p = 0; p < P; p++)
			{
				//look at the first layer
				cout << "\n";
				for (int i = 0; i < N; i++)
				{
					cout << " " << trainingValues[p][i];
				}
			}
		}

		B = new int[N];

		int * minn = new int[N]; for (int i = 0; i < N; i++) minn[i] = 1000000;
		int * maxn = new int[N]; for (int i = 0; i < N; i++) maxn[i] = -1000000;

		for (int p = 0; p < P; p++)
		{
			for (int i = 0; i < N; i++)
			{
				if (display) cout << "\n tr:" << trainingValues[p][i] << ", min=" << minn[i] << ", max=" << maxn[i];
				minn[i] = min(minn[i], trainingValues[p][i]);
				maxn[i] = max(maxn[i], trainingValues[p][i]);
			}
		}
		if (display)
		{
			cout << "\n min: ";
			for (int i = 0; i < N; i++) cout << " " << minn[i];
			cout << "\n max: ";
			for (int i = 0; i < N; i++) cout << " " << maxn[i];
		}

		for (int i = 0; i < N; i++) B[i] = minn[i];

		if (display)
		{
			cout << "\n B: ";
			for (int i = 0; i < N; i++) cout << " " << B[i];
		}
		D = new int[N];

		for (int i = 0; i < N; i++) D[i] = maxn[i] - minn[i] + 1;
		D[N - 1] = maxn[N - 1];
		if (display)
		{
			cout << "\n D: ";
			for (int i = 0; i < N; i++) cout << " " << D[i];
		}
		delete[] minn;
		delete[] maxn;

		int NoPA; //number of attributes for this particular attribute

		appendFirst(D[0], 0, -1);

		for (int k = 1; k<N - 1; k++)
		{
			appendInter(P, D[k], k, -1);
		}
		int x0;

		//let's fill the structure

		for (int p = 0; p<P; p++)
		{
			//look at the first layer
			x0 = trainingValues[p][0] - B[0];
			int nextIndex;
			layer * L = first;

			if (L->fuzzy[x0][0] == 100)
			{
				//there's one road here already, let's go down on it
				nextIndex = L->index[x0][0];
			}
			else
			{
				//no road right here, let's make one
				L->index[x0][0] = L->contents;
				nextIndex = L->contents;
				RO = ceil(D[0]*RP);
				if (display) cout << "R0=" << RO;
				calculateAndSetFuzzy1D(L, x0, 0, RO);
				L->fuzzy[x0][0] = 100;
				L->contents++;
			}
			//look at the rest of the layers


			for (int ii = 1; ii < N - 1; ii++)
			{
				L = L->next;
				x0 = trainingValues[p][ii] - B[ii];
				if (L->fuzzy[nextIndex][x0] == 100)
				{
					//there's one road here already, let's go down on it
					nextIndex = L->index[nextIndex][x0];
					if (!L->next) myfile >> ff; //redundant class, let's just read it and do nothing with it
				}
				else
				{
					//no road right here, let's make one
					if (L->next)
					{
						L->index[nextIndex][x0] = L->contents;
					}
					else
					{
						myfile >> ff;
						int y0 = trainingValues[p][N - 1] - B[N - 1];
						L->index[nextIndex][x0] = y0;
					}
					RO = ceil(D[ii] * RP);
					if (display) cout << "R0=" << RO;
					calculateAndSetFuzzy1D(L, nextIndex, x0, RO);

					L->fuzzy[nextIndex][x0] = 100;

					nextIndex = L->contents;
					L->contents++;
				}
			}
			//we reached and processed the last layer
		}

		end = clock(); 	elapsed_sec = double(end - start) / CLOCKS_PER_SEC; 	cout << "\n time spent on training: " << elapsed_sec;

		for (int i = 0; i<P; i++)
		{
			delete[] trainingValues[i];
		}
		delete[] trainingValues;

		cout << "\n Training is done!";
	}
	void addNumber2List(int id0)
	{		
		if (!numbersFirst)
		{
			numberNode * numbersFirst = new numberNode;
			numbersFirst->next = NULL;
			numbersFirst->data = id0;
			numbersLast = numbersFirst;
		}
		else
		{
			numberNode * t = new numberNode;
			numbersLast->next = t;
			t->next = NULL;
			t->data = id0;
			numbersLast = t;
		}
		numberCounter++;
	}
	void add2IDlist(IDnode * q, short id0)
	{
		IDnode * p;
		for (p = q; p->next; p = p->next) {}
		IDnode * t = new IDnode;
		p->next = t;
		t->next = NULL;
		t->ID = id0;
	}
	void appendList2List(IDnode * a, IDnode * b) //copy contents of b into a
	{
		IDnode * c = b;
		while (c)
		{
			add2IDlist(a, c->ID);
			c = c->next;
		}
	}
	void gatherItemsFromLayer(IDnode * IDN, layer * L, int *X, int index0, short paramNo)
	{
		int I0 = paramNo;//L->paramNo0;
		int nextIndex;
		if (L->next)
		{
			for (int i = 0; i < X[I0]; i++)
			{
				if (L->fuzzy[index0][i] > 0)
				{
					//there's one road here already, let's go down on it
					nextIndex = L->index[index0][i];
					gatherItemsFromLayer(IDN, L->next, X, nextIndex, paramNo + 1);
				}
			}
		}
		else
		{
			for (int i = 0; i < X[I0]; i++)
			{
				if (L->fuzzy[index0][i] > 0)
				{
					//there's one road here already, let's go down on it
					nextIndex = L->index[index0][i];
					appendList2List(IDN, items[nextIndex]);
				}
			}
		}
	}

	IDnode * gatherItems(short N, int *X)
	{
		IDnode * results = new IDnode;
		results->ID = -1;
		results->next = NULL;
		layer * L = first;
		int I0 = L->paramNo0;
		int nextIndex;
		for (int i = 0; i < X[0]; i++)
		{
			if (L->fuzzy[i][0] > 0)
			{
				//there's one road here already, let's go down on it
				nextIndex = L->index[i][0];
				gatherItemsFromLayer(results, L->next, X, nextIndex, 1);
			}
			else
			{
				//no road right here, nor in the vicinity, return with DUNNO
				//cout<<"\n returned at the first one";
				//return results;
			}
		}
		return results;
	}

	response evaluateCurrentSample(int N)
	{
		cout << "\n";
		response resp;
		resp.classID = -1;
		resp.fuzzyValue = 1;
		resp.probability = 0;
		resp.s = "";
		layer * L = first;
		int I0 = L->paramNo0;
		int I1;
		short x0 = currentSample[I0] * A[I0] - B[I0];
		short x1;
		int nextIndex;

		if (L->fuzzy[x0][0]>0)
		{
			//there's one road here already, let's go down on it
			nextIndex = L->index[x0][0];
			resp.fuzzyValue /= 100;
			cout <<"\n     "<< nextIndex << " @(" << 0 << "/" << x0 << ")";
		}
		else
		{
			//no road right here, nor in the vicinity, return with DUNNO
			cout << "\n     nothing @(" << 0 << "/" << x0 << ")";
			return resp;
		}

		//look at the rest of the layers
		L = L->next;
		bool l = true;
		while (L)
		{
			I0 = L->paramNo0;
			x0 = currentSample[I0]; 
			if (L->fuzzy[nextIndex][x0]>0)
			{
				//there's one road here already, let's go down on it
				if (L->index[nextIndex][x0] >= 0)
				{
					cout << "\n     " << L->index[nextIndex][x0] << " @(" << nextIndex << "/" << x0 << ")";
					
					if (!L->next)
					{
						resp.classID = L->index[nextIndex][x0];
						resp.s += " found class at";
						resp.s += L->ID;
						return resp;
					}
					nextIndex = L->index[nextIndex][x0];
				}
				else
				{
					resp.classID = L->index[nextIndex][x0];
					resp.s += " found class at";
					resp.s += L->ID;
					return resp;
				}
			}
			else
			{
				//no road right here, nor in the vicinity, return with what we know so far
				cout << ">\n     nothing @(" << I0 << "/" << x0 << ")";
				resp.s += "returned at ";
				resp.s += L->ID;
				return resp;
			}
			L = L->next;
		}

		//we reached the last, the so-called "class" layer
		resp.s += " went full length";
		return resp;
	}
	response evaluateCurrentSampleNoDisplay(int N)
	{
		response resp;
		resp.classID = -1;
		resp.fuzzyValue = 1;
		resp.probability = 0;
		resp.s = "";
		layer * L = first;
		int I0 = L->paramNo0;
		int I1;
		int x0 = round(currentSample[I0] * A[I0]) - B[I0];
		int x1;
		int nextIndex;

		if (L->fuzzy[x0][0]>0)
		{
			//there's one road here already, let's go down on it
			nextIndex = L->index[x0][0];
			resp.fuzzyValue /= 100;
		}
		else
		{
			//no road right here, nor in the vicinity, return with DUNNO
			return resp;
		}

		//look at the rest of the layers
		L = L->next;
		bool l = true;
		while (L)
		{
			I0 = L->paramNo0;
			x0 = round(currentSample[I0]*A[I0])-B[I0];
			if (L->fuzzy[nextIndex][x0]>0)
			{
				//there's one road here already, let's go down on it
				if (L->index[nextIndex][x0] >= 0)
				{
					if (!L->next)
					{
						resp.classID = L->index[nextIndex][x0];
						resp.s += " found class at";
						resp.s += L->ID;
						return resp;
					}
					nextIndex = L->index[nextIndex][x0];
				}
				else
				{
					resp.classID = L->index[nextIndex][x0];
					resp.s += " found class at";
					resp.s += L->ID;
					return resp;
				}
			}
			else
			{
				//no road right here, nor in the vicinity, return with what we know so far
				resp.s += "returned at ";
				resp.s += L->ID;
				return resp;
			}
			L = L->next;
		}

		//we reached the last, the so-called "class" layer
		resp.s += " went full length";
		return resp;
	}

	void evaluate(char * file)
	{
		int correct = 0;

		int N; int P;
		fstream myfile;
		myfile.open(file, ios::in);
	//	fstream g("log.dat", ios::out);
		myfile >> P;
		myfile >> N;
		currentSample = new double[N];

		//// Create the dynamic arrays for the temporary storage of the training data
		//double ** samples = new  double*[P];
		//for (int i=0; i<P; i++)
		//{
		//	samples[i] = new  double[N];
		//	for (int j=0; j<N; j++)
		//	{
		//		samples[i][j]=0;
		//	}
		//}
		//short * classes = new short[P];
		//for (int j=0; j<P; j++)
		//{
		//	classes[j]=0;
		//}

		// read the data from file
		double t; short t1;
		for (int p = 0; p<P; p++)
		{
			for (int j = 0; j<N - 1; j++)
			{
				myfile >> t;
				currentSample[j] = t;
			}
			myfile >> t1;

			response resp = evaluateCurrentSample(N);
			//cout<<"\n    resp="<<(short)resp.classID;
			cout << "\n" << p;
			if ((short)resp.classID == t1)
			{
				cout << ": CORRECT";
				correct++;
			}
			else
			{
				cout << ": MISTAKE, getting " << (short)resp.classID << " instead of " << t1 << ", " << resp.s.c_str();
			}
			//if ((short)resp.classID==t1) { correct++;} 
		}
		myfile.close();
	//	g.close();
		cout << "\n\n ------- CORRECT: " << correct << "/" << P << "  (" << 100 * correct / P << "%)";
	}
	void evaluateNoDisplay(char * file)
	{
		int correct = 0;

		clock_t start, end;
		double elapsed_sec;
		start = clock();

		int N; int P;
		fstream myfile;
		myfile.open(file, ios::in);
		myfile >> P;
		myfile >> N;
		currentSample = new double[N];
		double t; double t1;
		for (int j = 0; j<N; j++)
		{
			myfile >> t;
		}
		// read the data from file

		for (int p = 0; p<P; p++)
		{
			for (int j = 0; j<N - 1; j++)
			{
				myfile >> t;
				currentSample[j] = t;
			}
			myfile >> t1;

			response resp = evaluateCurrentSampleNoDisplay(N);
			if ((short)resp.classID == round(t1*A[N-1])-B[N-1])
			{
				correct++;
			}
			else
			{
			}
		}
		myfile.close();
		end = clock(); 	elapsed_sec = double(end - start) / CLOCKS_PER_SEC; 	cout << "\n time spent on evaluation: " << elapsed_sec;
		cout << "\n\n ------- CORRECT: " << correct << "/" << P << "  (" << 100 * correct / P << "%)";
	}



	~SFIT()
	{
		layer *q;
		if (first == NULL)
			return;

		while (first != NULL)
		{
			q = first->next;

			if (first)
			{
				//cout << "\n deleting "<<(int)first->ID;
				if (first->fuzzy) { for (int i = 0; i<first->sizeX; i++) delete[] first->fuzzy[i]; delete[] first->fuzzy; }

				if (first->index) { for (int i = 0; i<first->sizeX; i++) delete[] first->index[i]; delete[] first->index; }
				delete first;
				first = q;
			}

		}
		if (numbersFirst)
		{
			while (numbersFirst)
			{
				numberNode * t = numbersFirst->next;
				delete numbersFirst;
				numbersFirst = t;
			}
		}
		if (currentSample) delete currentSample;
		if (A) delete A;
		if (B) delete B;
	}

};
#endif