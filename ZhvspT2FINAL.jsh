// this creates graphs of nu vs. delta pt^2 for the three different elements


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


// this creates histograms of the zh data for the 3 elements' heavy and deuterium targets

for (int i = 0; i < 6; i++) {	
	String hname = "hpT2_" + i;
	histyArray[i] = new H1F().create(hname, zh_nBins, (DataVector)vecArray[i].get(0),zh_lower,zh_higher);
	
}


// this next part creates the equal statistical bins of nu (1st bin has the data with the lowest values of zh, 6th bin has the data with the highest values of zh)


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




TreeHipo[] treeNumber = {tree1, tree1, tree2, tree2, tree3, tree3};

int[] treeTarget = {0, 1, 0, 1, 0, 1};

int[] entriesArray = {entries1, entries1, entries2, entries2, entries3, entries3};



for (int i = 0; i < 6; i++) {
	treeNumber[i].reset();
	vecArray[6*i+6] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh<"+finalBinNum1Array[i],entriesArray[i]);
	
	treeNumber[i].reset();
	vecArray[6*i+7] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum1Array[i]+"&&zh<"+finalBinNum2Array[i],entriesArray[i]);

	treeNumber[i].reset();
	vecArray[6*i+8] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum2Array[i]+"&&zh<"+finalBinNum3Array[i],entriesArray[i]);

	treeNumber[i].reset();
	vecArray[6*i+9] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum3Array[i]+"&&zh<"+finalBinNum4Array[i],entriesArray[i]);

	treeNumber[i].reset();
	vecArray[6*i+10] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum4Array[i]+"&&zh<"+finalBinNum5Array[i],entriesArray[i]);

	treeNumber[i].reset();
	vecArray[6*i+11] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum5Array[i],entriesArray[i]);
}


int pT2_nBins = 100;
double pT2_lower = 0.0;
double pT2_higher = 2.0;


for (int i = 6; i < 42; i++) {	
	String hname = "hpT2_" + i;
	histyArray[i] = new H1F().create(hname, pT2_nBins, (DataVector)vecArray[i].get(0),pT2_lower,pT2_higher);
	
}


// determines the average value of pT^2 from each interval of zh
double[] histyMean = new double[42];
for (int i = 0; i < 42; i++) {   
	histyMean[i] = histyArray[i].getMean();
}



// this part determines the values of delta pT^2

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


// this part determines the error bars for the graphs
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



// the bins stuff is copied again because the code worked better this way

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




//*****okay onto the middle stuff******

// this finds the average value of zh within each bin to plot against the average value of delta pT^2

// middle 1

double[] targetBinMiddle1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	targetBinMiddle1Array[i] = histyArray[i].integral(0,finalBin1Array[i])/2;	
}



double[] highIntegralMiddle1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegralMiddle1Array[i] = 0.0;	
}

int[] highBinMiddle1Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBinMiddle1Array[i] = -1;	
}


int[] lowBinMiddle1Array = new int[6];


double[] highCheckMiddle1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheckMiddle1Array[i] = 0.0;	
}

double[] lowCheckMiddle1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheckMiddle1Array[i] = 0.0;	
}

int[] finalBinMiddle1Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBinMiddle1Array[i] = 0;	
}

double[] finalBinNumMiddle1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNumMiddle1Array[i] = 0.0;	
}



for (int i = 0; i < 6; i++) {

	while (highIntegralMiddle1Array[i] < targetBinMiddle1Array[i]) {
		highBinMiddle1Array[i]++;
		highIntegralMiddle1Array[i] = histyArray[i].integral(0,highBinMiddle1Array[i]);
	}

}



for (int i = 0; i < 6; i++) {	
	lowBinMiddle1Array[i] = highBinMiddle1Array[i] - 1;	
}



double[] lowIntegralMiddle1Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegralMiddle1Array[i] = histyArray[i].integral(0,lowBinMiddle1Array[i]);	
}


for (int i = 0; i < 6; i++) {
	highCheckMiddle1Array[i] = Math.abs(targetBinMiddle1Array[i] - highIntegralMiddle1Array[i]);


	lowCheckMiddle1Array[i] = Math.abs(targetBinMiddle1Array[i] - lowIntegralMiddle1Array[i]);



	if (highCheckMiddle1Array[i] < lowCheckMiddle1Array[i]) {
	finalBinNumMiddle1Array[i] = histyArray[i].getDataX(highBinMiddle1Array[i]);
	finalBinMiddle1Array[i] = highBinMiddle1Array[i];
	} else {
	finalBinNumMiddle1Array[i] = histyArray[i].getDataX(lowBinMiddle1Array[i]);
	finalBinMiddle1Array[i] = lowBinMiddle1Array[i];
	}
}


// start of Middle2


double[] targetBinMiddle2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	targetBinMiddle2Array[i] = histyArray[i].integral(startNum2Array[i],finalBin2Array[i])/2;	
}



double[] highIntegralMiddle2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegralMiddle2Array[i] = 0;	
}


int[] highBinMiddle2Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBinMiddle2Array[i] = finalBinMiddle1Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegralMiddle2Array[i] < targetBinMiddle2Array[i]) {
		highBinMiddle2Array[i]++;
		highIntegralMiddle2Array[i] = histyArray[i].integral(startNum2Array[i],highBinMiddle2Array[i]);
	}
}

int[] lowBinMiddle2Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBinMiddle2Array[i] = highBinMiddle2Array[i]-1;	
}



double[] lowIntegralMiddle2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegralMiddle2Array[i] = histyArray[i].integral(startNum2Array[i],lowBinMiddle2Array[i]);	
}

double[] highCheckMiddle2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheckMiddle2Array[i] = 0;	
}

double[] lowCheckMiddle2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheckMiddle2Array[i] = 0;	
}


int[] finalBinMiddle2Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBinMiddle2Array[i] = 0;	
}

double[] finalBinNumMiddle2Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNumMiddle2Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheckMiddle2Array[i] = Math.abs(targetBinMiddle2Array[i] - highIntegralMiddle2Array[i]);


	lowCheckMiddle2Array[i] = Math.abs(targetBinMiddle2Array[i] - lowIntegralMiddle2Array[i]);



	if (highCheckMiddle2Array[i] < lowCheckMiddle2Array[i]) {
	finalBinNumMiddle2Array[i] = histyArray[i].getDataX(highBinMiddle2Array[i]);
	finalBinMiddle2Array[i] = highBinMiddle2Array[i];
	} else {
	finalBinNumMiddle2Array[i] = histyArray[i].getDataX(lowBinMiddle2Array[i]);
	finalBinMiddle2Array[i] = lowBinMiddle2Array[i];
	}
}


// start of Middle 3

double[] targetBinMiddle3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	targetBinMiddle3Array[i] = histyArray[i].integral(startNum3Array[i],finalBin3Array[i])/2;	
}



double[] highIntegralMiddle3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegralMiddle3Array[i] = 0;	
}


int[] highBinMiddle3Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBinMiddle3Array[i] = finalBinMiddle2Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegralMiddle3Array[i] < targetBinMiddle3Array[i]) {
		highBinMiddle3Array[i]++;
		highIntegralMiddle3Array[i] = histyArray[i].integral(startNum3Array[i],highBinMiddle3Array[i]);
	}
}

int[] lowBinMiddle3Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBinMiddle3Array[i] = highBinMiddle3Array[i]-1;	
}



double[] lowIntegralMiddle3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegralMiddle3Array[i] = histyArray[i].integral(startNum3Array[i],lowBinMiddle3Array[i]);	
}

double[] highCheckMiddle3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheckMiddle3Array[i] = 0;	
}

double[] lowCheckMiddle3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheckMiddle3Array[i] = 0;	
}


int[] finalBinMiddle3Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBinMiddle3Array[i] = 0;	
}

double[] finalBinNumMiddle3Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNumMiddle3Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheckMiddle3Array[i] = Math.abs(targetBinMiddle3Array[i] - highIntegralMiddle3Array[i]);


	lowCheckMiddle3Array[i] = Math.abs(targetBinMiddle3Array[i] - lowIntegralMiddle3Array[i]);



	if (highCheckMiddle3Array[i] < lowCheckMiddle3Array[i]) {
	finalBinNumMiddle3Array[i] = histyArray[i].getDataX(highBinMiddle3Array[i]);
	finalBinMiddle3Array[i] = highBinMiddle3Array[i];
	} else {
	finalBinNumMiddle3Array[i] = histyArray[i].getDataX(lowBinMiddle3Array[i]);
	finalBinMiddle3Array[i] = lowBinMiddle3Array[i];
	}
}

// now onto 4

double[] targetBinMiddle4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	targetBinMiddle4Array[i] = histyArray[i].integral(startNum4Array[i],finalBin4Array[i])/2;	
}



double[] highIntegralMiddle4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegralMiddle4Array[i] = 0;	
}


int[] highBinMiddle4Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBinMiddle4Array[i] = finalBinMiddle3Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegralMiddle4Array[i] < targetBinMiddle4Array[i]) {
		highBinMiddle4Array[i]++;
		highIntegralMiddle4Array[i] = histyArray[i].integral(startNum4Array[i],highBinMiddle4Array[i]);
	}
}



int[] lowBinMiddle4Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBinMiddle4Array[i] = highBinMiddle4Array[i]-1;	
}



double[] lowIntegralMiddle4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegralMiddle4Array[i] = histyArray[i].integral(startNum4Array[i],lowBinMiddle4Array[i]);	
}

double[] highCheckMiddle4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheckMiddle4Array[i] = 0;	
}

double[] lowCheckMiddle4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheckMiddle4Array[i] = 0;	
}


int[] finalBinMiddle4Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBinMiddle4Array[i] = 0;	
}

double[] finalBinNumMiddle4Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNumMiddle4Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheckMiddle4Array[i] = Math.abs(targetBinMiddle4Array[i] - highIntegralMiddle4Array[i]);


	lowCheckMiddle4Array[i] = Math.abs(targetBinMiddle4Array[i] - lowIntegralMiddle4Array[i]);



	if (highCheckMiddle4Array[i] < lowCheckMiddle4Array[i]) {
	finalBinNumMiddle4Array[i] = histyArray[i].getDataX(highBinMiddle4Array[i]);
	finalBinMiddle4Array[i] = highBinMiddle4Array[i];
	} else {
	finalBinNumMiddle4Array[i] = histyArray[i].getDataX(lowBinMiddle4Array[i]);
	finalBinMiddle4Array[i] = lowBinMiddle4Array[i];
	}
}



// now onto 5

double[] targetBinMiddle5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	targetBinMiddle5Array[i] = histyArray[i].integral(startNum5Array[i],finalBin5Array[i])/2;	
}



double[] highIntegralMiddle5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegralMiddle5Array[i] = 0;	
}


int[] highBinMiddle5Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBinMiddle5Array[i] = finalBinMiddle4Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegralMiddle5Array[i] < targetBinMiddle5Array[i]) {
		highBinMiddle5Array[i]++;
		highIntegralMiddle5Array[i] = histyArray[i].integral(startNum5Array[i],highBinMiddle5Array[i]);
	}
}



int[] lowBinMiddle5Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBinMiddle5Array[i] = highBinMiddle5Array[i]-1;	
}



double[] lowIntegralMiddle5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegralMiddle5Array[i] = histyArray[i].integral(startNum5Array[i],lowBinMiddle5Array[i]);	
}

double[] highCheckMiddle5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheckMiddle5Array[i] = 0;	
}

double[] lowCheckMiddle5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheckMiddle5Array[i] = 0;	
}


int[] finalBinMiddle5Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBinMiddle5Array[i] = 0;	
}

double[] finalBinNumMiddle5Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNumMiddle5Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheckMiddle5Array[i] = Math.abs(targetBinMiddle5Array[i] - highIntegralMiddle5Array[i]);


	lowCheckMiddle5Array[i] = Math.abs(targetBinMiddle5Array[i] - lowIntegralMiddle5Array[i]);



	if (highCheckMiddle5Array[i] < lowCheckMiddle5Array[i]) {
	finalBinNumMiddle5Array[i] = histyArray[i].getDataX(highBinMiddle5Array[i]);
	finalBinMiddle5Array[i] = highBinMiddle5Array[i];
	} else {
	finalBinNumMiddle5Array[i] = histyArray[i].getDataX(lowBinMiddle5Array[i]);
	finalBinMiddle5Array[i] = lowBinMiddle5Array[i];
	}
}

// now onto middle 6

int[] startNum6Array = new int[6];
for (int i = 0; i < 6; i++) {	
	startNum6Array[i] = finalBin5Array[i]+1;	
}

int[] finalBin6Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBin6Array[i] = 99;
}


double[] targetBinMiddle6Array = new double[6];
for (int i = 0; i < 6; i++) {	
	targetBinMiddle6Array[i] = histyArray[i].integral(startNum6Array[i],finalBin6Array[i])/2;	
}



double[] highIntegralMiddle6Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highIntegralMiddle6Array[i] = 0;	
}


int[] highBinMiddle6Array = new int[6];
for (int i = 0; i < 6; i++) {	
	highBinMiddle6Array[i] = finalBinMiddle5Array[i];	
}

for (int i = 0; i < 6; i++) {
	while (highIntegralMiddle6Array[i] < targetBinMiddle6Array[i]) {
		highBinMiddle6Array[i]++;
		highIntegralMiddle6Array[i] = histyArray[i].integral(startNum6Array[i],highBinMiddle6Array[i]);
	}
}



int[] lowBinMiddle6Array = new int[6];
for (int i = 0; i < 6; i++) {	
	lowBinMiddle6Array[i] = highBinMiddle6Array[i]-1;	
}



double[] lowIntegralMiddle6Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowIntegralMiddle6Array[i] = histyArray[i].integral(startNum6Array[i],lowBinMiddle6Array[i]);	
}

double[] highCheckMiddle6Array = new double[6];
for (int i = 0; i < 6; i++) {	
	highCheckMiddle6Array[i] = 0;	
}

double[] lowCheckMiddle6Array = new double[6];
for (int i = 0; i < 6; i++) {	
	lowCheckMiddle6Array[i] = 0;	
}


int[] finalBinMiddle6Array = new int[6];
for (int i = 0; i < 6; i++) {	
	finalBinMiddle6Array[i] = 0;	
}

double[] finalBinNumMiddle6Array = new double[6];
for (int i = 0; i < 6; i++) {	
	finalBinNumMiddle6Array[i] = 0;	
}


for (int i = 0; i < 6; i++) {
	highCheckMiddle6Array[i] = Math.abs(targetBinMiddle6Array[i] - highIntegralMiddle6Array[i]);


	lowCheckMiddle6Array[i] = Math.abs(targetBinMiddle6Array[i] - lowIntegralMiddle6Array[i]);



	if (highCheckMiddle6Array[i] < lowCheckMiddle6Array[i]) {
	finalBinNumMiddle6Array[i] = histyArray[i].getDataX(highBinMiddle6Array[i]);
	finalBinMiddle6Array[i] = highBinMiddle6Array[i];
	} else {
	finalBinNumMiddle6Array[i] = histyArray[i].getDataX(lowBinMiddle6Array[i]);
	finalBinMiddle6Array[i] = lowBinMiddle6Array[i];
	}
}



// this part is because there is a different average value of zh for each element bin: the deuterium average and the heavy target average


double[] middleNumberBin1Array = new double[3];
for (int i = 0; i < 3; i++) {	
	middleNumberBin1Array[i] = (finalBinNumMiddle1Array[2*i]+finalBinNumMiddle1Array[2*i+1])/2;	
}


double[] middleNumberBin2Array = new double[3];
for (int i = 0; i < 3; i++) {	
	middleNumberBin2Array[i] = (finalBinNumMiddle2Array[2*i]+finalBinNumMiddle2Array[2*i+1])/2;	
}

double[] middleNumberBin3Array = new double[3];
for (int i = 0; i < 3; i++) {	
	middleNumberBin3Array[i] = (finalBinNumMiddle3Array[2*i]+finalBinNumMiddle3Array[2*i+1])/2;	
}

double[] middleNumberBin4Array = new double[3];
for (int i = 0; i < 3; i++) {	
	middleNumberBin4Array[i] = (finalBinNumMiddle4Array[2*i]+finalBinNumMiddle4Array[2*i+1])/2;	
}


double[] middleNumberBin5Array = new double[3];
for (int i = 0; i < 3; i++) {	
	middleNumberBin5Array[i] = (finalBinNumMiddle5Array[2*i]+finalBinNumMiddle5Array[2*i+1])/2;	
}

double[] middleNumberBin6Array = new double[3];
for (int i = 0; i < 3; i++) {	
	middleNumberBin6Array[i] = (finalBinNumMiddle6Array[2*i]+finalBinNumMiddle6Array[2*i+1])/2;	
}

// end of middle stuff



GraphErrors graph = new GraphErrors();
graph.addPoint(middleNumberBin1Array[0],carbon_differencebin1,0,deltacarbondiffmeanbin1);
graph.addPoint(middleNumberBin2Array[0],carbon_differencebin2,0,deltacarbondiffmeanbin2);
graph.addPoint(middleNumberBin3Array[0],carbon_differencebin3,0,deltacarbondiffmeanbin3);
graph.addPoint(middleNumberBin4Array[0],carbon_differencebin4,0,deltacarbondiffmeanbin4);
graph.addPoint(middleNumberBin5Array[0],carbon_differencebin5,0,deltacarbondiffmeanbin5);
graph.addPoint(middleNumberBin6Array[0],carbon_differencebin6,0,deltacarbondiffmeanbin6);
graph.setTitle("#Delta pT^2 vs zh for Carbon");
graph.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph.setTitleX("zh");


GraphErrors graph1 = new GraphErrors();
graph1.addPoint(middleNumberBin1Array[1],iron_differencebin1,0,deltairondiffmeanbin1);
graph1.addPoint(middleNumberBin2Array[1],iron_differencebin2,0,deltairondiffmeanbin2);
graph1.addPoint(middleNumberBin3Array[1],iron_differencebin3,0,deltairondiffmeanbin3);
graph1.addPoint(middleNumberBin4Array[1],iron_differencebin4,0,deltairondiffmeanbin4);
graph1.addPoint(middleNumberBin5Array[1],iron_differencebin5,0,deltairondiffmeanbin5);
graph1.addPoint(middleNumberBin6Array[1],iron_differencebin6,0,deltairondiffmeanbin6);
graph1.setTitle("#Delta pT^2 vs zh for Iron");
graph1.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph1.setTitleX("zh");


GraphErrors graph2 = new GraphErrors();
graph2.addPoint(middleNumberBin1Array[2],lead_differencebin1,0,deltaleaddiffmeanbin1);
graph2.addPoint(middleNumberBin2Array[2],lead_differencebin2,0,deltaleaddiffmeanbin2);
graph2.addPoint(middleNumberBin3Array[2],lead_differencebin3,0,deltaleaddiffmeanbin3);
graph2.addPoint(middleNumberBin4Array[2],lead_differencebin4,0,deltaleaddiffmeanbin4);
graph2.addPoint(middleNumberBin5Array[2],lead_differencebin5,0,deltaleaddiffmeanbin5);
graph2.addPoint(middleNumberBin6Array[2],lead_differencebin6,0,deltaleaddiffmeanbin6);
graph2.setTitle("#Delta pT^2 vs zh for Lead");
graph2.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph2.setTitleX("zh");



// I'm still not sure what this part with the directory does

String dirname = "/electron";
TDirectory dir = new TDirectory();
dir.mkdir(dirname);
dir.cd(dirname);



int c1a_title_size = 30;
TCanvas c1 = new TCanvas("c1",800,800);
// number columns number rows
c1.divide(3,1); 


c1.draw(graph);


c1.cd(1);
c1.draw(graph1);


c1.cd(2);
c1.draw(graph2);


