// this is the 6 graphs code from ZhCounts8.2 with the general code improvements from ZhCounts 6.3 -- basically ZhCounts8.2 with code fixed

// this is basiclally the same as 8.2 except I created variables for the graph parameters (like pT_lower, etc)

// to do: I could maybe create more variables for the graph points (like the 12, etc)
 
// I created algoritihm code stuff for getting the differences for the elements in 9.2 but it's not that important to do


// commented out the ugly vec stuff

import org.jlab.jnp.hipo4.data.*;
import org.jlab.jnp.hipo4.data.*;
import org.jlab.jnp.hipo4.io.*;
import org.jlab.clas.physics.*;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.*;
import org.jlab.groot.ui.*;
import org.jlab.groot.tree.*; 
import javax.swing.JFrame;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.math.*;
import java.lang.Math;

int maxEvents = -1;
String fileName1 = "ntuple_C_NposCOPY.hipo"; 
String fileName2 = "ntuple_Fe_NposCOPY.hipo"; 
String fileName3 = "ntuple_Pb_NposCOPY.hipo"; 



TreeHipo tree1 = new TreeHipo(fileName1,"protonTree::tree"); 
int entries1 = tree1.getEntries();   
System.out.println(" ENTRIES = " + entries1);  
TreeHipo tree2 = new TreeHipo(fileName2,"protonTree::tree"); 
int entries2 = tree2.getEntries();   
System.out.println(" ENTRIES = " + entries2);  
TreeHipo tree3 = new TreeHipo(fileName3,"protonTree::tree"); 
int entries3 = tree3.getEntries();   
System.out.println(" ENTRIES = " + entries3);  



List vec = tree1.getDataVectors("zh","pFidCut==1&&eFidCut==1&&iTgt==0",entries1);
List vec1 = tree1.getDataVectors("zh","pFidCut==1&&eFidCut==1&&iTgt==1",entries1);
List vec2 = tree2.getDataVectors("zh","pFidCut==1&&eFidCut==1&&iTgt==0",entries2);
List vec3 = tree2.getDataVectors("zh","pFidCut==1&&eFidCut==1&&iTgt==1",entries2);
List vec4 = tree3.getDataVectors("zh","pFidCut==1&&eFidCut==1&&iTgt==0",entries3);
List vec5 = tree3.getDataVectors("zh","pFidCut==1&&eFidCut==1&&iTgt==1",entries3);



H1F[] histyArray = new H1F[42];

List[] vecArray = new List[42];
vecArray[0] = vec;
vecArray[1] = vec1;
vecArray[2] = vec2;
vecArray[3] = vec3;
vecArray[4] = vec4;
vecArray[5] = vec5;



int zh_nBins = 100;
double zh_lower = 0.0;
double zh_higher = 1.25;

for (int i = 0; i < 6; i++) {	
	String hname = "hpT2_" + i;
	histyArray[i] = new H1F().create(hname, zh_nBins, (DataVector)vecArray[i].get(0),zh_lower,zh_higher);
	
}


double[] targetBinArray = new double[6];
for (int i = 0; i < 6; i++) {	
	targetBinArray[i] = histyArray[i].getIntegral()/6;	
}

double[] highIntegral1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegral1Array[i] = 0.0;	
}

int[] highBin1Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBin1Array[i] = -1;	
}


int[] lowBin1Array = new int[6];


double[] highCheck1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheck1Array[i] = 0.0;	
}

double[] lowCheck1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheck1Array[i] = 0.0;	
}

int[] finalBin1Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBin1Array[i] = 0;	
}

double[] finalBinNum1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNum1Array[i] = 0.0;	
}


for (int i = 0; i < 6; i++) {

	while (highIntegral1Array[i] < targetBinArray[i]) {
		highBin1Array[i]++;
		highIntegral1Array[i] = histyArray[i].integral(0,highBin1Array[i]);
	}

}



for (int i = 0; i < 6; i++) {	
	lowBin1Array[i] = highBin1Array[i] - 1;	
}



double[] lowIntegral1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegral1Array[i] = histyArray[i].integral(0,lowBin1Array[i]);	
}


for (int i = 0; i < 6; i++) {
	highCheck1Array[i] = Math.abs(targetBinArray[i] - highIntegral1Array[i]);


	lowCheck1Array[i] = Math.abs(targetBinArray[i] - lowIntegral1Array[i]);



	if (highCheck1Array[i] < lowCheck1Array[i]) {
	finalBinNum1Array[i] = histyArray[i].getDataX(highBin1Array[i]);
	finalBin1Array[i] = highBin1Array[i];
	} else {
	finalBinNum1Array[i] = histyArray[i].getDataX(lowBin1Array[i]);
	finalBin1Array[i] = lowBin1Array[i];
	}
}

// number 2

int[] startNum2Array = new int[6];
for (int i = 0; i < 6; i++) {	
	startNum2Array[i] = finalBin1Array[i]+1;	
}

double[] highIntegral2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegral2Array[i] = 0;	
}


int[] highBin2Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBin2Array[i] = finalBin1Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegral2Array[i] < targetBinArray[i]) {
		highBin2Array[i]++;
		highIntegral2Array[i] = histyArray[i].integral(startNum2Array[i],highBin2Array[i]);
	}
}

int[] lowBin2Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBin2Array[i] = highBin2Array[i]-1;	
}



double[] lowIntegral2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegral2Array[i] = histyArray[i].integral(startNum2Array[i],lowBin2Array[i]);	
}

double[] highCheck2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheck2Array[i] = 0;	
}

double[] lowCheck2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheck2Array[i] = 0;	
}


int[] finalBin2Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBin2Array[i] = 0;	
}

double[] finalBinNum2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNum2Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheck2Array[i] = Math.abs(targetBinArray[i] - highIntegral2Array[i]);


	lowCheck2Array[i] = Math.abs(targetBinArray[i] - lowIntegral2Array[i]);



	if (highCheck2Array[i] < lowCheck2Array[i]) {
	finalBinNum2Array[i] = histyArray[i].getDataX(highBin2Array[i]);
	finalBin2Array[i] = highBin2Array[i];
	} else {
	finalBinNum2Array[i] = histyArray[i].getDataX(lowBin2Array[i]);
	finalBin2Array[i] = lowBin2Array[i];
	}
}

// now onto number 3

int[] startNum3Array = new int[6];
for (int i = 0; i < 6; i++) {	
	startNum3Array[i] = finalBin2Array[i]+1;	
}

double[] highIntegral3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegral3Array[i] = 0;	
}


int[] highBin3Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBin3Array[i] = finalBin2Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegral3Array[i] < targetBinArray[i]) {
		highBin3Array[i]++;
		highIntegral3Array[i] = histyArray[i].integral(startNum3Array[i],highBin3Array[i]);
	}
}

int[] lowBin3Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBin3Array[i] = highBin3Array[i]-1;	
}



double[] lowIntegral3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegral3Array[i] = histyArray[i].integral(startNum3Array[i],lowBin3Array[i]);	
}

double[] highCheck3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheck3Array[i] = 0;	
}

double[] lowCheck3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheck3Array[i] = 0;	
}


int[] finalBin3Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBin3Array[i] = 0;	
}

double[] finalBinNum3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNum3Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheck3Array[i] = Math.abs(targetBinArray[i] - highIntegral3Array[i]);


	lowCheck3Array[i] = Math.abs(targetBinArray[i] - lowIntegral3Array[i]);



	if (highCheck3Array[i] < lowCheck3Array[i]) {
	finalBinNum3Array[i] = histyArray[i].getDataX(highBin3Array[i]);
	finalBin3Array[i] = highBin3Array[i];
	} else {
	finalBinNum3Array[i] = histyArray[i].getDataX(lowBin3Array[i]);
	finalBin3Array[i] = lowBin3Array[i];
	}
}

// now onto number 4

int[] startNum4Array = new int[6];
for (int i = 0; i < 6; i++) {	
	startNum4Array[i] = finalBin3Array[i]+1;	
}

double[] highIntegral4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegral4Array[i] = 0;	
}


int[] highBin4Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBin4Array[i] = finalBin3Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegral4Array[i] < targetBinArray[i]) {
		highBin4Array[i]++;
		highIntegral4Array[i] = histyArray[i].integral(startNum4Array[i],highBin4Array[i]);
	}
}

int[] lowBin4Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBin4Array[i] = highBin4Array[i]-1;	
}



double[] lowIntegral4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegral4Array[i] = histyArray[i].integral(startNum4Array[i],lowBin4Array[i]);	
}

double[] highCheck4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheck4Array[i] = 0;	
}

double[] lowCheck4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheck4Array[i] = 0;	
}


int[] finalBin4Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBin4Array[i] = 0;	
}

double[] finalBinNum4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNum4Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheck4Array[i] = Math.abs(targetBinArray[i] - highIntegral4Array[i]);


	lowCheck4Array[i] = Math.abs(targetBinArray[i] - lowIntegral4Array[i]);



	if (highCheck4Array[i] < lowCheck4Array[i]) {
	finalBinNum4Array[i] = histyArray[i].getDataX(highBin4Array[i]);
	finalBin4Array[i] = highBin4Array[i];
	} else {
	finalBinNum4Array[i] = histyArray[i].getDataX(lowBin4Array[i]);
	finalBin4Array[i] = lowBin4Array[i];
	}
}

// now onto 5


int[] startNum5Array = new int[6];
for (int i = 0; i < 6; i++) {	
	startNum5Array[i] = finalBin4Array[i]+1;	
}

double[] highIntegral5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegral5Array[i] = 0;	
}


int[] highBin5Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBin5Array[i] = finalBin4Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegral5Array[i] < targetBinArray[i]) {
		highBin5Array[i]++;
		highIntegral5Array[i] = histyArray[i].integral(startNum5Array[i],highBin5Array[i]);
	}
}

int[] lowBin5Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBin5Array[i] = highBin5Array[i]-1;	
}



double[] lowIntegral5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegral5Array[i] = histyArray[i].integral(startNum5Array[i],lowBin5Array[i]);	
}

double[] highCheck5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheck5Array[i] = 0;	
}

double[] lowCheck5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheck5Array[i] = 0;	
}


int[] finalBin5Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBin5Array[i] = 0;	
}

double[] finalBinNum5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNum5Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheck5Array[i] = Math.abs(targetBinArray[i] - highIntegral5Array[i]);


	lowCheck5Array[i] = Math.abs(targetBinArray[i] - lowIntegral5Array[i]);



	if (highCheck5Array[i] < lowCheck5Array[i]) {
	finalBinNum5Array[i] = histyArray[i].getDataX(highBin5Array[i]);
	finalBin5Array[i] = highBin5Array[i];
	} else {
	finalBinNum5Array[i] = histyArray[i].getDataX(lowBin5Array[i]);
	finalBin5Array[i] = lowBin5Array[i];
	}
}



// 

TreeHipo[] treeNumber = {tree1, tree1, tree2, tree2, tree3, tree3};

int[] treeTarget = {0, 1, 0, 1, 0, 1};

int[] entriesArray = {entries1, entries1, entries2, entries2, entries3, entries3};

// !!!!!!!!!!!!!!! make this into an array?


/*
vecArray[6] = vec6;
vecArray[7] = vec7;
vecArray[8] = vec8;
vecArray[9] = vec9;
vecArray[10] = vec10;
vecArray[11] = vec11;


vecArray[12] = vec12;
vecArray[13] = vec13;
vecArray[14] = vec14;
vecArray[15] = vec15;
vecArray[16] = vec16;
vecArray[17] = vec17;


vecArray[18] = vec18;
vecArray[19] = vec19;
vecArray[20] = vec20;
vecArray[21] = vec21;
vecArray[22] = vec22;
vecArray[23] = vec23;


vecArray[24] = vec24;
vecArray[25] = vec25;
vecArray[26] = vec26;
vecArray[27] = vec27;
vecArray[28] = vec28;
vecArray[29] = vec29;

vecArray[30] = vec30;
vecArray[31] = vec31;
vecArray[32] = vec32;
vecArray[33] = vec33;
vecArray[34] = vec34;
vecArray[35] = vec35;

vecArray[36] = vec36;
vecArray[37] = vec37;
vecArray[38] = vec38;
vecArray[39] = vec39;
vecArray[40] = vec40;
vecArray[41] = vec41;

*/

for (int i = 0; i < 6; i++) {
	treeNumber[i].reset();
	vecArray[6*i+6] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh<"+finalBinNum1Array[i],entriesArray[i]);
	
	treeNumber[i].reset();
	vecArray[6*i+7] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum1Array[i]+"&&zh<"+finalBinNum2Array[i],entriesArray[i]);

	treeNumber[i].reset();
	vecArray[6*i+8] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum2Array[i]+"&&zh<"+finalBinNum3Array[i],entriesArray[i]);

	treeNumber[i].reset();
	vecArray[6*i+9] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum3Array[i]+"&&zh<"+finalBinNum4Array[i],entriesArray[i]);

// 
	treeNumber[i].reset();
	vecArray[6*i+10] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum4Array[i]+"&&zh<"+finalBinNum5Array[i],entriesArray[i]);

	treeNumber[i].reset();
	vecArray[6*i+11] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum5Array[i],entriesArray[i]);
}


int pT2_nBins = 100;
double pT2_lower = 0.0;
double pT2_higher = 2.0;


for (int i = 0; i < 42; i++) {	
	String hname = "hpT2_" + i;
	histyArray[i] = new H1F().create(hname, pT2_nBins, (DataVector)vecArray[i].get(0),pT2_lower,pT2_higher);
	
}



double[] histyMean = new double[42];
for (int i = 0; i < 42; i++) {   
	histyMean[i] = histyArray[i].getMean();
	//System.out.println("number " + i + " " + histyMean[i]);
}


/*
double[] histyIntegral = new double[42];
for (int i = 0; i < 42; i++) {   
	histyIntegral[i] = histyArray[i].getIntegral();
	System.out.println("number " + i + " " + histyIntegral[i]);
}

*/

// meaning of deut higher than heavy target?


double carbon_differencebin1 = histyMean[6]-histyMean[12]; 
double iron_differencebin1 = histyMean[18] - histyMean[24]; 
double lead_differencebin1 = histyMean[30] - histyMean[36]; 


double carbon_differencebin2 = histyMean[7]-histyMean[13]; 
double iron_differencebin2 = histyMean[19] - histyMean[25]; 
double lead_differencebin2 = histyMean[31] - histyMean[37]; 

double carbon_differencebin3 = histyMean[8]-histyMean[14]; 
double iron_differencebin3 = histyMean[20] - histyMean[26]; 
double lead_differencebin3 = histyMean[32] - histyMean[38]; 


double carbon_differencebin4 = histyMean[9]-histyMean[15]; 
double iron_differencebin4 = histyMean[21] - histyMean[27]; 
double lead_differencebin4 = histyMean[33] - histyMean[39]; 


double carbon_differencebin5 = histyMean[16]-histyMean[10]; 
double iron_differencebin5 = histyMean[22] - histyMean[28]; 
double lead_differencebin5 = histyMean[34] - histyMean[40]; 


double carbon_differencebin6 = histyMean[17]-histyMean[11]; 
double iron_differencebin6 = histyMean[29] - histyMean[23]; 
double lead_differencebin6 = histyMean[41] - histyMean[35]; 


double[] errorMean = new double[42];
for (int i = 0; i < 42; i++) {	
	double N_d = 0; 
	double err_d = 0;
	double counts_d = 0;
	double xi_d = 0;
	double dN2_d = 0;
	double dN_d = 0;
	double CT_d = 0;
	double dCT_d = 0;
	double e1_d = 0;
	double dmean_d = 0;
	for(int j=0;j<100;j++){
		err_d = histyArray[i].getBinError(j);
		counts_d = histyArray[i].getBinContent(j);
		xi_d = histyArray[i].getDataX(j);
		N_d = N_d + xi_d*counts_d;
		dN2_d = xi_d*xi_d*err_d*err_d;
		}
	dN_d = Math.sqrt(dN2_d);
	CT_d = histyArray[i].getEntries();
	dCT_d = Math.sqrt(CT_d);
	e1_d = (dN_d/N_d)*(dN_d/N_d)+(dCT_d/CT_d)*(dCT_d/CT_d);
	dmean_d = histyMean[i]*Math.sqrt(e1_d);
	errorMean[i] = dmean_d;
}




double deltacarbondiffmeanbin1 = 0;
double deltairondiffmeanbin1 = 0;
double deltaleaddiffmeanbin1 = 0;


deltacarbondiffmeanbin1 = Math.sqrt((errorMean[6])*(errorMean[6])-(errorMean[12])*(errorMean[12]));
deltairondiffmeanbin1 = Math.sqrt((errorMean[18])*(errorMean[18])-(errorMean[24])*(errorMean[24]));
deltaleaddiffmeanbin1 = Math.sqrt((errorMean[30])*(errorMean[30])-(errorMean[36])*(errorMean[36])); 


double deltacarbondiffmeanbin2 = 0;
double deltairondiffmeanbin2 = 0;
double deltaleaddiffmeanbin2 = 0;


deltacarbondiffmeanbin2 = Math.sqrt((errorMean[7])*(errorMean[7])-(errorMean[13])*(errorMean[13]));
deltairondiffmeanbin2 = Math.sqrt((errorMean[19])*(errorMean[19])-(errorMean[25])*(errorMean[25]));
deltaleaddiffmeanbin2 = Math.sqrt((errorMean[31])*(errorMean[31])-(errorMean[37])*(errorMean[37])); 


double deltacarbondiffmeanbin3 = 0;
double deltairondiffmeanbin3 = 0;
double deltaleaddiffmeanbin3 = 0;

deltacarbondiffmeanbin3 = Math.sqrt((errorMean[8])*(errorMean[8])-(errorMean[14])*(errorMean[14]));
deltairondiffmeanbin3 = Math.sqrt((errorMean[20])*(errorMean[20])-(errorMean[26])*(errorMean[26]));
deltaleaddiffmeanbin3 = Math.sqrt((errorMean[32])*(errorMean[32])-(errorMean[38])*(errorMean[38])); 


double deltacarbondiffmeanbin4 = 0;
double deltairondiffmeanbin4 = 0;
double deltaleaddiffmeanbin4 = 0;

deltacarbondiffmeanbin4 = Math.sqrt((errorMean[9])*(errorMean[9])-(errorMean[15])*(errorMean[15]));
deltairondiffmeanbin4 = Math.sqrt((errorMean[21])*(errorMean[21])-(errorMean[27])*(errorMean[27]));
deltaleaddiffmeanbin4 = Math.sqrt((errorMean[33])*(errorMean[33])-(errorMean[39])*(errorMean[39])); 


double deltacarbondiffmeanbin5 = 0;
double deltairondiffmeanbin5 = 0;
double deltaleaddiffmeanbin5 = 0;

deltacarbondiffmeanbin5 = Math.sqrt((errorMean[10])*(errorMean[10])-(errorMean[16])*(errorMean[16]));
deltairondiffmeanbin5 = Math.sqrt((errorMean[22])*(errorMean[22])-(errorMean[28])*(errorMean[28]));
deltaleaddiffmeanbin5 = Math.sqrt((errorMean[34])*(errorMean[34])-(errorMean[40])*(errorMean[40])); 



double deltacarbondiffmeanbin6 = 0;
double deltairondiffmeanbin6 = 0;
double deltaleaddiffmeanbin6 = 0;

deltacarbondiffmeanbin6 = Math.sqrt((errorMean[11])*(errorMean[11])-(errorMean[17])*(errorMean[17]));
deltairondiffmeanbin6 = Math.sqrt((errorMean[23])*(errorMean[23])-(errorMean[29])*(errorMean[29]));
deltaleaddiffmeanbin6 = Math.sqrt((errorMean[35])*(errorMean[35])-(errorMean[41])*(errorMean[41])); 


GraphErrors graph = new GraphErrors();
graph.addPoint(12,carbon_differencebin2,0,deltacarbondiffmeanbin2);
graph.addPoint(56,iron_differencebin2,0,deltairondiffmeanbin2);
graph.addPoint(208,lead_differencebin2,0,deltaleaddiffmeanbin2);
graph.setTitle("#Delta pT^2 vs Atomic Number for second interval");
graph.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph.setTitleX("Atomic Number");


GraphErrors graph1 = new GraphErrors();
graph1.addPoint(12,carbon_differencebin3,0,deltacarbondiffmeanbin3);
graph1.addPoint(56,iron_differencebin3,0,deltairondiffmeanbin3);
graph1.addPoint(208,lead_differencebin3,0,deltaleaddiffmeanbin3);
graph1.setTitle("#Delta pT^2 vs Atomic Number for third interval");
graph1.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph1.setTitleX("Atomic Number");


GraphErrors graph2 = new GraphErrors();
graph2.addPoint(12,carbon_differencebin4,0,deltacarbondiffmeanbin4);
graph2.addPoint(56,iron_differencebin4,0,deltairondiffmeanbin4);
graph2.addPoint(208,lead_differencebin4,0,deltaleaddiffmeanbin4);
graph2.setTitle("#Delta pT^2 vs Atomic Number for fourth interval");
graph2.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph2.setTitleX("Atomic Number");

GraphErrors graph3 = new GraphErrors();
graph3.addPoint(12,carbon_differencebin1,0,deltacarbondiffmeanbin1);
graph3.addPoint(56,iron_differencebin1,0,deltairondiffmeanbin1);
graph3.addPoint(208,lead_differencebin1,0,deltaleaddiffmeanbin1);
graph3.setTitle("#Delta pT^2 vs Atomic Number for first interval");
graph3.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph3.setTitleX("Atomic Number");


GraphErrors graph4 = new GraphErrors();
graph4.addPoint(12,carbon_differencebin5,0,deltacarbondiffmeanbin5);
graph4.addPoint(56,iron_differencebin5,0,deltairondiffmeanbin5);
graph4.addPoint(208,lead_differencebin5,0,deltaleaddiffmeanbin5);
graph4.setTitle("#Delta pT^2 vs Atomic Number for fifth interval");
graph4.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph4.setTitleX("Atomic Number");


GraphErrors graph5 = new GraphErrors();
graph5.addPoint(12,carbon_differencebin6,0,deltacarbondiffmeanbin6);
graph5.addPoint(56,iron_differencebin6,0,deltairondiffmeanbin6);
graph5.addPoint(208,lead_differencebin6,0,deltaleaddiffmeanbin6);
graph5.setTitle("#Delta pT^2 vs Atomic Number for sixth interval");
graph5.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph5.setTitleX("Atomic Number");



// I'm still not sure what this part with the directory does

String dirname = "/electron";
TDirectory dir = new TDirectory();
dir.mkdir(dirname);
dir.cd(dirname);



//F1D func = new F1D("func","[a]+[b]*x", 0, 210);
//DataFitter.fit(func, graph, "Q");
F1D func1 = new F1D("func1","[a]+[b]*x+[c]*x*x", 0, 210);
DataFitter.fit(func1, graph, "Q");
//func.setLineColor(6);
//func.setLineWidth(3);
func1.setLineColor(7);
func1.setLineWidth(3);

//F1D func2 = new F1D("func2","[a]+[b]*x", 0, 210);
//DataFitter.fit(func2, graph1, "Q");
F1D func3 = new F1D("func3","[a]+[b]*x+[c]*x*x", 0, 210);
DataFitter.fit(func3, graph1, "Q");
//func2.setLineColor(6);
//func2.setLineWidth(3);
func3.setLineColor(7);
func3.setLineWidth(3);

//F1D func4 = new F1D("func4","[a]+[b]*x", 0, 210);
//DataFitter.fit(func4, graph2, "Q");
F1D func5 = new F1D("func5","[a]+[b]*x+[c]*x*x", 0, 210);
DataFitter.fit(func5, graph2, "Q");
//func4.setLineColor(6);
//func4.setLineWidth(3);
func5.setLineColor(7);
func5.setLineWidth(3);

//F1D func6 = new F1D("func6","[a]+[b]*x", 0, 210);
//DataFitter.fit(func6, graph3, "Q");
F1D func7 = new F1D("func7","[a]+[b]*x+[c]*x*x", 0, 210);
DataFitter.fit(func7, graph3, "Q");
//func6.setLineColor(6);
//func6.setLineWidth(3);
func7.setLineColor(7);
func7.setLineWidth(3);


F1D func9 = new F1D("func9","[a]+[b]*x+[c]*x*x", 0, 210);
DataFitter.fit(func9, graph4, "Q");
func9.setLineColor(7);
func9.setLineWidth(3);


F1D func11 = new F1D("func11","[a]+[b]*x+[c]*x*x", 0, 210);
DataFitter.fit(func11, graph5, "Q");
func11.setLineColor(7);
func11.setLineWidth(3);


int c1a_title_size = 30;
TCanvas c1 = new TCanvas("c1",800,800);
// number columns number rows
c1.divide(3,2); 



c1.draw(graph3);
c1.draw(func7, "same");

c1.cd(1);
//c1.draw(histyArray[2]);
//dir.addDataSet(histyArray[2]);
c1.draw(graph);
c1.draw(func1, "same");

c1.cd(2);
//c1.draw(histyArray[2]);
//dir.addDataSet(histyArray[2]);
c1.draw(graph1);
c1.draw(func3, "same");


c1.cd(3);
//c1.draw(histyArray[2]);
//dir.addDataSet(histyArray[2]);
c1.draw(graph2);
c1.draw(func5, "same");


c1.cd(4);
//c1.draw(histyArray[2]);
//dir.addDataSet(histyArray[2]);
c1.draw(graph4);
c1.draw(func9, "same");


c1.cd(5);
//c1.draw(histyArray[2]);
//dir.addDataSet(histyArray[2]);
c1.draw(graph5);
c1.draw(func11, "same");
