// creates 6 graphs of delta pT^2 vs. A divided by bins of zh. On each graph is the data further divided into 3 bins of q2

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

// ntuple files to use
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



H1F[] histyArray = new H1F[50];

List[] vecArray = new List[50];
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


// this next part creates 6 equal statistical bins of zh (1st bin has the data with the lowest values of zh, 6th bin has the data with the highest values of zh)
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

List[] pT2vecArray = new List[36];

for (int i = 0; i < 6; i++) {
	treeNumber[i].reset();
	pT2vecArray[6*i] = treeNumber[i].getDataVectors("q2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh<"+finalBinNum1Array[i],entriesArray[i]);
	
	treeNumber[i].reset();
	pT2vecArray[6*i+1] = treeNumber[i].getDataVectors("q2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum1Array[i]+"&&zh<"+finalBinNum2Array[i],entriesArray[i]);

	treeNumber[i].reset();
	pT2vecArray[6*i+2] = treeNumber[i].getDataVectors("q2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum2Array[i]+"&&zh<"+finalBinNum3Array[i],entriesArray[i]);

	treeNumber[i].reset();
	pT2vecArray[6*i+3] = treeNumber[i].getDataVectors("q2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum3Array[i]+"&&zh<"+finalBinNum4Array[i],entriesArray[i]);

// 
	treeNumber[i].reset();
	pT2vecArray[6*i+4] = treeNumber[i].getDataVectors("q2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum4Array[i]+"&&zh<"+finalBinNum5Array[i],entriesArray[i]);

	treeNumber[i].reset();
	pT2vecArray[6*i+5] = treeNumber[i].getDataVectors("q2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&zh>"+finalBinNum5Array[i],entriesArray[i]);
}



int q2_nBins = 100;
double q2_lower = 1.0;
double q2_higher = 4.1;

H1F[] pT2histyArray = new H1F[36];




for (int i = 0; i < 36; i++) {	
	String h1stname = "h1stpT2_" + i;
	pT2histyArray[i] = new H1F().create(h1stname, q2_nBins, (DataVector)pT2vecArray[i].get(0),q2_lower,q2_higher);
	
}



// this next part looks at each individual histogram created ^^^^ there (36 histograms are 6 bins of zh for each of the heavy elements and deuterium elements) it then divides the q2 data in each histogram into 3 bins


double[] a1sttargetBinArray = new double[36];
for (int i = 0; i < 36; i++) {	
	a1sttargetBinArray[i] = pT2histyArray[i].getIntegral()/3;	
}

double[] a1sthighIntegral1Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1sthighIntegral1Array[i] = 0.0;	
}

int[] a1sthighBin1Array = new int[36];
for (int i = 0; i < 36; i++) {	
	a1sthighBin1Array[i] = -1;	
}


int[] a1stlowBin1Array = new int[36];


double[] a1sthighCheck1Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1sthighCheck1Array[i] = 0.0;	
}

double[] a1stlowCheck1Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1stlowCheck1Array[i] = 0.0;	
}

int[] a1stfinalBin1Array = new int[36];
for (int i = 0; i < 36; i++) {	
	a1stfinalBin1Array[i] = 0;	
}

double[] a1stfinalBinNum1Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1stfinalBinNum1Array[i] = 0.0;	
}


for (int i = 0; i < 36; i++) {

	while (a1sthighIntegral1Array[i] < a1sttargetBinArray[i]) {
		a1sthighBin1Array[i]++;
		a1sthighIntegral1Array[i] = pT2histyArray[i].integral(0,a1sthighBin1Array[i]);
	}

}



for (int i = 0; i < 36; i++) {	
	a1stlowBin1Array[i] = a1sthighBin1Array[i] - 1;	
}



double[] a1stlowIntegral1Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1stlowIntegral1Array[i] = pT2histyArray[i].integral(0,a1stlowBin1Array[i]);	
}


for (int i = 0; i < 36; i++) {
	a1sthighCheck1Array[i] = Math.abs(a1sttargetBinArray[i] - a1sthighIntegral1Array[i]);


	a1stlowCheck1Array[i] = Math.abs(a1sttargetBinArray[i] - a1stlowIntegral1Array[i]);



	if (a1sthighCheck1Array[i] < a1stlowCheck1Array[i]) {
	a1stfinalBinNum1Array[i] = pT2histyArray[i].getDataX(a1sthighBin1Array[i]);
	a1stfinalBin1Array[i] = a1sthighBin1Array[i];
	} else {
	a1stfinalBinNum1Array[i] = pT2histyArray[i].getDataX(a1stlowBin1Array[i]);
	a1stfinalBin1Array[i] = a1stlowBin1Array[i];
	}
}

// number 2

int[] a1ststartNum2Array = new int[36];
for (int i = 0; i < 36; i++) {	
	a1ststartNum2Array[i] = a1stfinalBin1Array[i]+1;	
}

double[] a1sthighIntegral2Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1sthighIntegral2Array[i] = 0;	
}


int[] a1sthighBin2Array = new int[36];
for (int i = 0; i < 36; i++) {	
	a1sthighBin2Array[i] = a1stfinalBin1Array[i];	
}

for (int i = 0; i < 36; i++) {
	while (a1sthighIntegral2Array[i] < a1sttargetBinArray[i]) {
		a1sthighBin2Array[i]++;
		a1sthighIntegral2Array[i] = pT2histyArray[i].integral(a1ststartNum2Array[i],a1sthighBin2Array[i]);
	}
}

int[] a1stlowBin2Array = new int[36];
for (int i = 0; i < 36; i++) {	
	a1stlowBin2Array[i] = a1sthighBin2Array[i]-1;	
}



double[] a1stlowIntegral2Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1stlowIntegral2Array[i] = pT2histyArray[i].integral(a1ststartNum2Array[i],a1stlowBin2Array[i]);	
}

double[] a1sthighCheck2Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1sthighCheck2Array[i] = 0;	
}

double[] a1stlowCheck2Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1stlowCheck2Array[i] = 0;	
}


int[] a1stfinalBin2Array = new int[36];
for (int i = 0; i < 36; i++) {	
	a1stfinalBin2Array[i] = 0;	
}

double[] a1stfinalBinNum2Array = new double[36];
for (int i = 0; i < 36; i++) {	
	a1stfinalBinNum2Array[i] = 0;	
}


for (int i = 0; i < 36; i++) {
	a1sthighCheck2Array[i] = Math.abs(a1sttargetBinArray[i] - a1sthighIntegral2Array[i]);


	a1stlowCheck2Array[i] = Math.abs(a1sttargetBinArray[i] - a1stlowIntegral2Array[i]);



	if (a1sthighCheck2Array[i] < a1stlowCheck2Array[i]) {
	a1stfinalBinNum2Array[i] = pT2histyArray[i].getDataX(a1sthighBin2Array[i]);
	a1stfinalBin2Array[i] = a1sthighBin2Array[i];
	} else {
	a1stfinalBinNum2Array[i] = pT2histyArray[i].getDataX(a1stlowBin2Array[i]);
	a1stfinalBin2Array[i] = a1stlowBin2Array[i];
	}
}




List[] finalvecArray = new List[108];



for (int i = 0; i < 6; i++) {
	
// these are the 3 q2s for zh 1
	treeNumber[i].reset();
	finalvecArray[18*i] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2<"+a1stfinalBinNum1Array[6*i]+"&&zh<"+finalBinNum1Array[i],entriesArray[i]);
	

	treeNumber[i].reset();
	finalvecArray[18*i+1] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum1Array[6*i]+"&&q2<"+a1stfinalBinNum2Array[6*i]+"&&zh<"+finalBinNum1Array[i],entriesArray[i]);


	treeNumber[i].reset();
	finalvecArray[18*i+2] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum2Array[6*i]+"&&zh<"+finalBinNum1Array[i],entriesArray[i]);

// these are the 3 q2s for zh 2

	treeNumber[i].reset();
	finalvecArray[18*i+3] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2<"+a1stfinalBinNum1Array[6*i+1]+"&&zh>"+finalBinNum1Array[i]+"&&zh<"+finalBinNum2Array[i],entriesArray[i]);
	

	treeNumber[i].reset();
	finalvecArray[18*i+4] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum1Array[6*i+1]+"&&q2<"+a1stfinalBinNum2Array[6*i+1]+"&&zh>"+finalBinNum1Array[i]+"&&zh<"+finalBinNum2Array[i],entriesArray[i]);


	treeNumber[i].reset();
	finalvecArray[18*i+5] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum2Array[6*i+1]+"&&zh>"+finalBinNum1Array[i]+"&&zh<"+finalBinNum2Array[i],entriesArray[i]);


// these are the 3 q2s for zh 3

	treeNumber[i].reset();
	finalvecArray[18*i+6] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2<"+a1stfinalBinNum1Array[6*i+2]+"&&zh>"+finalBinNum2Array[i]+"&&zh<"+finalBinNum3Array[i],entriesArray[i]);
	

	treeNumber[i].reset();
	finalvecArray[18*i+7] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum1Array[6*i+2]+"&&q2<"+a1stfinalBinNum2Array[6*i+2]+"&&zh>"+finalBinNum2Array[i]+"&&zh<"+finalBinNum3Array[i],entriesArray[i]);


	treeNumber[i].reset();
	finalvecArray[18*i+8] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum2Array[6*i+2]+"&&zh>"+finalBinNum2Array[i]+"&&zh<"+finalBinNum3Array[i],entriesArray[i]);


// these are the 3 q2s for zh 4

	treeNumber[i].reset();
	finalvecArray[18*i+9] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2<"+a1stfinalBinNum1Array[6*i+3]+"&&zh>"+finalBinNum3Array[i]+"&&zh<"+finalBinNum4Array[i],entriesArray[i]);
	

	treeNumber[i].reset();
	finalvecArray[18*i+10] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum1Array[6*i+3]+"&&q2<"+a1stfinalBinNum2Array[6*i+3]+"&&zh>"+finalBinNum3Array[i]+"&&zh<"+finalBinNum4Array[i],entriesArray[i]);


	treeNumber[i].reset();
	finalvecArray[18*i+11] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum2Array[6*i+3]+"&&zh>"+finalBinNum3Array[i]+"&&zh<"+finalBinNum4Array[i],entriesArray[i]);


// these are the 3 q2s for zh 5

	treeNumber[i].reset();
	finalvecArray[18*i+12] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2<"+a1stfinalBinNum1Array[6*i+4]+"&&zh>"+finalBinNum4Array[i]+"&&zh<"+finalBinNum5Array[i],entriesArray[i]);
	

	treeNumber[i].reset();
	finalvecArray[18*i+13] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum1Array[6*i+4]+"&&q2<"+a1stfinalBinNum2Array[6*i+4]+"&&zh>"+finalBinNum4Array[i]+"&&zh<"+finalBinNum5Array[i],entriesArray[i]);


	treeNumber[i].reset();
	finalvecArray[18*i+14] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum2Array[6*i+4]+"&&zh>"+finalBinNum4Array[i]+"&&zh<"+finalBinNum5Array[i],entriesArray[i]);

// these are the 3 q2s for zh 6

	treeNumber[i].reset();
	finalvecArray[18*i+15] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2<"+a1stfinalBinNum1Array[6*i+5]+"&&zh>"+finalBinNum5Array[i],entriesArray[i]);
	

	treeNumber[i].reset();
	finalvecArray[18*i+16] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum1Array[6*i+5]+"&&q2<"+a1stfinalBinNum2Array[6*i+5]+"&&zh>"+finalBinNum5Array[i],entriesArray[i]);


	treeNumber[i].reset();
	finalvecArray[18*i+17] = treeNumber[i].getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt=="+treeTarget[i]+"&&q2>"+a1stfinalBinNum2Array[6*i+5]+"&&zh>"+finalBinNum5Array[i],entriesArray[i]);


}



int pT2_nBins = 100;
double pT2_lower = 0.0;
double pT2_higher = 2.0;


H1F[] finalhistyArray = new H1F[108];


for (int i = 0; i < 108; i++) {	
	String hname = "hpT2_" + i;
	finalhistyArray[i] = new H1F().create(hname, pT2_nBins, (DataVector)finalvecArray[i].get(0),pT2_lower,pT2_higher);
	
}



// determines the average value of pT^2 from each zh&q2 division
double[] histyMean = new double[108];
for (int i = 0; i < 108; i++) {   
	histyMean[i] = finalhistyArray[i].getMean();
	
}


// this part determines the values of delta pT^2
double carbon_differencebin1q1 = histyMean[0]-histyMean[18]; 
double carbon_differencebin1q2 = histyMean[1]-histyMean[19]; 
double carbon_differencebin1q3 = histyMean[2]-histyMean[20]; 


double iron_differencebin1q1 = histyMean[36] - histyMean[54]; 
double iron_differencebin1q2 = histyMean[37] - histyMean[55]; 
double iron_differencebin1q3 = histyMean[38] - histyMean[56]; 


double lead_differencebin1q1 = histyMean[72] - histyMean[90]; 
double lead_differencebin1q2 = histyMean[73] - histyMean[91]; 
double lead_differencebin1q3 = histyMean[74] - histyMean[92]; 



double carbon_differencebin2q1 = histyMean[3]-histyMean[21]; 
double carbon_differencebin2q2 = histyMean[4]-histyMean[22]; 
double carbon_differencebin2q3 = histyMean[5]-histyMean[23]; 


double iron_differencebin2q1 = histyMean[39] - histyMean[57]; 
double iron_differencebin2q2 = histyMean[40] - histyMean[58]; 
double iron_differencebin2q3 = histyMean[41] - histyMean[59]; 


double lead_differencebin2q1 = histyMean[75] - histyMean[93]; 
double lead_differencebin2q2 = histyMean[76] - histyMean[94]; 
double lead_differencebin2q3 = histyMean[77] - histyMean[95]; 


double carbon_differencebin3q1 = histyMean[6]-histyMean[24]; 
double carbon_differencebin3q2 = histyMean[7]-histyMean[25]; 
double carbon_differencebin3q3 = histyMean[8]-histyMean[26]; 


double iron_differencebin3q1 = histyMean[42] - histyMean[60]; 
double iron_differencebin3q2 = histyMean[43] - histyMean[61]; 
double iron_differencebin3q3 = histyMean[44] - histyMean[62]; 


double lead_differencebin3q1 = histyMean[78] - histyMean[96]; 
double lead_differencebin3q2 = histyMean[79] - histyMean[97]; 
double lead_differencebin3q3 = histyMean[80] - histyMean[98]; 


double carbon_differencebin4q1 = histyMean[9]-histyMean[27]; 
double carbon_differencebin4q2 = histyMean[10]-histyMean[28]; 
double carbon_differencebin4q3 = histyMean[29]-histyMean[11]; 


double iron_differencebin4q1 = histyMean[45] - histyMean[63]; 
double iron_differencebin4q2 = histyMean[46] - histyMean[64]; 
double iron_differencebin4q3 = histyMean[47] - histyMean[65]; 


double lead_differencebin4q1 = histyMean[81] - histyMean[99]; 
double lead_differencebin4q2 = histyMean[82] - histyMean[100]; 
double lead_differencebin4q3 = histyMean[83] - histyMean[101]; 


double carbon_differencebin5q1 = histyMean[12]-histyMean[30]; 
double carbon_differencebin5q2 = histyMean[13]-histyMean[31]; 
double carbon_differencebin5q3 = histyMean[32]-histyMean[14];  


double iron_differencebin5q1 = histyMean[48] - histyMean[66]; 
double iron_differencebin5q2 = histyMean[49] - histyMean[67]; 
double iron_differencebin5q3 = histyMean[68] - histyMean[50]; 



double lead_differencebin5q1 = histyMean[84] - histyMean[102]; 
double lead_differencebin5q2 = histyMean[85] - histyMean[103]; 
double lead_differencebin5q3 = histyMean[86] - histyMean[104]; 



double carbon_differencebin6q1 = histyMean[33]-histyMean[15];
double carbon_differencebin6q2 = histyMean[34]-histyMean[16]; 
double carbon_differencebin6q3 = histyMean[35]-histyMean[17]; 

 
double iron_differencebin6q1 = histyMean[69] - histyMean[51]; 
double iron_differencebin6q2 = histyMean[70] - histyMean[52]; 
double iron_differencebin6q3 = histyMean[71] - histyMean[53]; 



double lead_differencebin6q1 = histyMean[105] - histyMean[87]; 
double lead_differencebin6q2 = histyMean[106] - histyMean[88]; 
double lead_differencebin6q3 = histyMean[107] - histyMean[89]; 


// this part determines the error bars for the graphs
double[] errorMean = new double[108];
for (int i = 0; i < 108; i++) {	
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
		err_d = finalhistyArray[i].getBinError(j);
		counts_d = finalhistyArray[i].getBinContent(j);
		xi_d = finalhistyArray[i].getDataX(j);
		N_d = N_d + xi_d*counts_d;
		dN2_d = xi_d*xi_d*err_d*err_d;
		}
	dN_d = Math.sqrt(dN2_d);
	CT_d = finalhistyArray[i].getEntries();
	dCT_d = Math.sqrt(CT_d);
	e1_d = (dN_d/N_d)*(dN_d/N_d)+(dCT_d/CT_d)*(dCT_d/CT_d);
	dmean_d = histyMean[i]*Math.sqrt(e1_d);
	errorMean[i] = dmean_d;
}



double deltacarbondiffmeanbin1q1 = 0;
double deltacarbondiffmeanbin1q2 = 0;
double deltacarbondiffmeanbin1q3 = 0;


double deltairondiffmeanbin1q1 = 0;
double deltairondiffmeanbin1q2 = 0;
double deltairondiffmeanbin1q3 = 0;


double deltaleaddiffmeanbin1q1 = 0;
double deltaleaddiffmeanbin1q2 = 0;
double deltaleaddiffmeanbin1q3 = 0;


deltacarbondiffmeanbin1q1 = Math.sqrt((errorMean[0])*(errorMean[0])-(errorMean[18])*(errorMean[18]));
deltacarbondiffmeanbin1q2 = Math.sqrt((errorMean[1])*(errorMean[1])-(errorMean[19])*(errorMean[19]));
deltacarbondiffmeanbin1q3 = Math.sqrt((errorMean[2])*(errorMean[2])-(errorMean[20])*(errorMean[20]));



deltairondiffmeanbin1q1 = Math.sqrt((errorMean[36])*(errorMean[36])-(errorMean[54])*(errorMean[54]));
deltairondiffmeanbin1q2 = Math.sqrt((errorMean[37])*(errorMean[37])-(errorMean[55])*(errorMean[55]));
deltairondiffmeanbin1q3 = Math.sqrt((errorMean[38])*(errorMean[38])-(errorMean[56])*(errorMean[56]));



deltaleaddiffmeanbin1q1 = Math.sqrt((errorMean[72])*(errorMean[72])-(errorMean[90])*(errorMean[90])); 
deltaleaddiffmeanbin1q2 = Math.sqrt((errorMean[73])*(errorMean[73])-(errorMean[91])*(errorMean[91])); 
deltaleaddiffmeanbin1q3 = Math.sqrt((errorMean[74])*(errorMean[74])-(errorMean[92])*(errorMean[92])); 


double deltacarbondiffmeanbin2q1 = 0;
double deltacarbondiffmeanbin2q2 = 0;
double deltacarbondiffmeanbin2q3 = 0;


double deltairondiffmeanbin2q1 = 0;
double deltairondiffmeanbin2q2 = 0;
double deltairondiffmeanbin2q3 = 0;



double deltaleaddiffmeanbin2q1 = 0;
double deltaleaddiffmeanbin2q2 = 0;
double deltaleaddiffmeanbin2q3 = 0;



deltacarbondiffmeanbin2q1 = Math.sqrt((errorMean[3])*(errorMean[3])-(errorMean[21])*(errorMean[21]));
deltacarbondiffmeanbin2q2 = Math.sqrt((errorMean[4])*(errorMean[4])-(errorMean[22])*(errorMean[22]));
deltacarbondiffmeanbin2q3 = Math.sqrt((errorMean[5])*(errorMean[5])-(errorMean[23])*(errorMean[23]));



deltairondiffmeanbin2q1 = Math.sqrt((errorMean[39])*(errorMean[39])-(errorMean[57])*(errorMean[57]));
deltairondiffmeanbin2q2 = Math.sqrt((errorMean[40])*(errorMean[40])-(errorMean[58])*(errorMean[58]));
deltairondiffmeanbin2q3 = Math.sqrt((errorMean[41])*(errorMean[41])-(errorMean[59])*(errorMean[59]));


deltaleaddiffmeanbin2q1 = Math.sqrt((errorMean[75])*(errorMean[75])-(errorMean[93])*(errorMean[93])); 
deltaleaddiffmeanbin2q2 = Math.sqrt((errorMean[76])*(errorMean[76])-(errorMean[94])*(errorMean[94])); 
deltaleaddiffmeanbin2q3 = Math.sqrt((errorMean[77])*(errorMean[77])-(errorMean[95])*(errorMean[95])); 



double deltacarbondiffmeanbin3q1 = 0;
double deltacarbondiffmeanbin3q2 = 0;
double deltacarbondiffmeanbin3q3 = 0;



double deltairondiffmeanbin3q1 = 0;
double deltairondiffmeanbin3q2 = 0;
double deltairondiffmeanbin3q3 = 0;



double deltaleaddiffmeanbin3q1 = 0;
double deltaleaddiffmeanbin3q2 = 0;
double deltaleaddiffmeanbin3q3 = 0;


deltacarbondiffmeanbin3q1 = Math.sqrt((errorMean[6])*(errorMean[6])-(errorMean[24])*(errorMean[24]));
deltacarbondiffmeanbin3q2 = Math.sqrt((errorMean[7])*(errorMean[7])-(errorMean[25])*(errorMean[25]));
deltacarbondiffmeanbin3q3 = Math.sqrt((errorMean[8])*(errorMean[8])-(errorMean[26])*(errorMean[26]));



deltairondiffmeanbin3q1 = Math.sqrt((errorMean[42])*(errorMean[42])-(errorMean[60])*(errorMean[60]));
deltairondiffmeanbin3q2 = Math.sqrt((errorMean[43])*(errorMean[43])-(errorMean[61])*(errorMean[61]));
deltairondiffmeanbin3q3 = Math.sqrt((errorMean[44])*(errorMean[44])-(errorMean[62])*(errorMean[62]));


deltaleaddiffmeanbin3q1 = Math.sqrt((errorMean[78])*(errorMean[78])-(errorMean[96])*(errorMean[96])); 
deltaleaddiffmeanbin3q2 = Math.sqrt((errorMean[79])*(errorMean[79])-(errorMean[97])*(errorMean[97]));
deltaleaddiffmeanbin3q3 = Math.sqrt((errorMean[80])*(errorMean[80])-(errorMean[98])*(errorMean[98]));



double deltacarbondiffmeanbin4q1 = 0;
double deltacarbondiffmeanbin4q2 = 0;
double deltacarbondiffmeanbin4q3 = 0;



double deltairondiffmeanbin4q1 = 0;
double deltairondiffmeanbin4q2 = 0;
double deltairondiffmeanbin4q3 = 0;



double deltaleaddiffmeanbin4q1 = 0;
double deltaleaddiffmeanbin4q2 = 0;
double deltaleaddiffmeanbin4q3 = 0;


deltacarbondiffmeanbin4q1 = Math.sqrt((errorMean[9])*(errorMean[9])-(errorMean[27])*(errorMean[27]));
deltacarbondiffmeanbin4q2 = Math.sqrt((errorMean[10])*(errorMean[10])-(errorMean[28])*(errorMean[28]));
deltacarbondiffmeanbin4q3 = Math.sqrt((errorMean[11])*(errorMean[11])-(errorMean[29])*(errorMean[29]));



deltairondiffmeanbin4q1 = Math.sqrt((errorMean[45])*(errorMean[45])-(errorMean[63])*(errorMean[63]));
deltairondiffmeanbin4q2 = Math.sqrt((errorMean[46])*(errorMean[46])-(errorMean[64])*(errorMean[64]));
deltairondiffmeanbin4q3 = Math.sqrt((errorMean[47])*(errorMean[47])-(errorMean[65])*(errorMean[65]));



deltaleaddiffmeanbin4q1 = Math.sqrt((errorMean[81])*(errorMean[81])-(errorMean[99])*(errorMean[99])); 
deltaleaddiffmeanbin4q2 = Math.sqrt((errorMean[82])*(errorMean[82])-(errorMean[100])*(errorMean[100])); 
deltaleaddiffmeanbin4q3 = Math.sqrt((errorMean[83])*(errorMean[83])-(errorMean[101])*(errorMean[101])); 



double deltacarbondiffmeanbin5q1 = 0;
double deltacarbondiffmeanbin5q2 = 0;
double deltacarbondiffmeanbin5q3 = 0;



double deltairondiffmeanbin5q1 = 0;
double deltairondiffmeanbin5q2 = 0;
double deltairondiffmeanbin5q3 = 0;



double deltaleaddiffmeanbin5q1 = 0;
double deltaleaddiffmeanbin5q2 = 0;
double deltaleaddiffmeanbin5q3 = 0;


deltacarbondiffmeanbin5q1 = Math.sqrt((errorMean[12])*(errorMean[12])-(errorMean[30])*(errorMean[30]));
deltacarbondiffmeanbin5q2 = Math.sqrt((errorMean[13])*(errorMean[13])-(errorMean[31])*(errorMean[31]));
deltacarbondiffmeanbin5q3 = Math.sqrt((errorMean[14])*(errorMean[14])-(errorMean[32])*(errorMean[32]));



deltairondiffmeanbin5q1 = Math.sqrt((errorMean[48])*(errorMean[48])-(errorMean[66])*(errorMean[66]));
deltairondiffmeanbin5q2 = Math.sqrt((errorMean[49])*(errorMean[49])-(errorMean[67])*(errorMean[67]));
deltairondiffmeanbin5q3 = Math.sqrt((errorMean[50])*(errorMean[50])-(errorMean[68])*(errorMean[68]));



deltaleaddiffmeanbin5q1 = Math.sqrt((errorMean[84])*(errorMean[84])-(errorMean[102])*(errorMean[102])); 
deltaleaddiffmeanbin5q2 = Math.sqrt((errorMean[85])*(errorMean[85])-(errorMean[103])*(errorMean[103])); 
deltaleaddiffmeanbin5q3 = Math.sqrt((errorMean[86])*(errorMean[86])-(errorMean[104])*(errorMean[104])); 



double deltacarbondiffmeanbin6q1 = 0;
double deltacarbondiffmeanbin6q2 = 0;
double deltacarbondiffmeanbin6q3 = 0;



double deltairondiffmeanbin6q1 = 0;
double deltairondiffmeanbin6q2 = 0;
double deltairondiffmeanbin6q3 = 0;


double deltaleaddiffmeanbin6q1 = 0;
double deltaleaddiffmeanbin6q2 = 0;
double deltaleaddiffmeanbin6q3 = 0;


deltacarbondiffmeanbin6q1 = Math.sqrt((errorMean[15])*(errorMean[15])-(errorMean[33])*(errorMean[33]));
deltacarbondiffmeanbin6q2 = Math.sqrt((errorMean[16])*(errorMean[16])-(errorMean[34])*(errorMean[34]));
deltacarbondiffmeanbin6q3 = Math.sqrt((errorMean[17])*(errorMean[17])-(errorMean[35])*(errorMean[35]));


deltairondiffmeanbin6q1 = Math.sqrt((errorMean[51])*(errorMean[51])-(errorMean[69])*(errorMean[69]));
deltairondiffmeanbin6q2 = Math.sqrt((errorMean[52])*(errorMean[52])-(errorMean[70])*(errorMean[70]));
deltairondiffmeanbin6q3 = Math.sqrt((errorMean[53])*(errorMean[53])-(errorMean[71])*(errorMean[71]));


deltaleaddiffmeanbin6q1 = Math.sqrt((errorMean[87])*(errorMean[87])-(errorMean[105])*(errorMean[105])); 
deltaleaddiffmeanbin6q2 = Math.sqrt((errorMean[88])*(errorMean[88])-(errorMean[106])*(errorMean[106])); 
deltaleaddiffmeanbin6q3 = Math.sqrt((errorMean[107])*(errorMean[107])-(errorMean[89])*(errorMean[89])); 





GraphErrors graphq1 = new GraphErrors();
GraphErrors graphq2 = new GraphErrors();
GraphErrors graphq3 = new GraphErrors();


graphq1.addPoint(12,carbon_differencebin2q1,0,deltacarbondiffmeanbin2q1);
graphq2.addPoint(12,carbon_differencebin2q2,0,deltacarbondiffmeanbin2q2);
graphq3.addPoint(12,carbon_differencebin2q3,0,deltacarbondiffmeanbin2q3);



graphq1.addPoint(56,iron_differencebin2q1,0,deltairondiffmeanbin2q1);
graphq2.addPoint(56,iron_differencebin2q2,0,deltairondiffmeanbin2q2);
graphq3.addPoint(56,iron_differencebin2q3,0,deltairondiffmeanbin2q3);


graphq1.addPoint(208,lead_differencebin2q1,0,deltaleaddiffmeanbin2q1);
graphq2.addPoint(208,lead_differencebin2q2,0,deltaleaddiffmeanbin2q2);
graphq3.addPoint(208,lead_differencebin2q3,0,deltaleaddiffmeanbin2q3);


graphq1.setTitle("#Delta pT^2 vs Atomic Number for 2nd interval of zh");
graphq1.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graphq1.setTitleX("Atomic Number");


GraphErrors graph1q1 = new GraphErrors();
GraphErrors graph1q2 = new GraphErrors();
GraphErrors graph1q3 = new GraphErrors();


graph1q1.addPoint(12,carbon_differencebin3q1,0,deltacarbondiffmeanbin3q1);
graph1q2.addPoint(12,carbon_differencebin3q2,0,deltacarbondiffmeanbin3q2);
graph1q3.addPoint(12,carbon_differencebin3q3,0,deltacarbondiffmeanbin3q3);


graph1q1.addPoint(56,iron_differencebin3q1,0,deltairondiffmeanbin3q1);
graph1q2.addPoint(56,iron_differencebin3q2,0,deltairondiffmeanbin3q2);
graph1q3.addPoint(56,iron_differencebin3q3,0,deltairondiffmeanbin3q3);


graph1q1.addPoint(208,lead_differencebin3q1,0,deltaleaddiffmeanbin3q1);
graph1q2.addPoint(208,lead_differencebin3q2,0,deltaleaddiffmeanbin3q2);
graph1q3.addPoint(208,lead_differencebin3q3,0,deltaleaddiffmeanbin3q3);


graph1q1.setTitle("#Delta pT^2 vs Atomic Number for 3rd interval of zh");
graph1q1.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph1q1.setTitleX("Atomic Number");


GraphErrors graph2q1 = new GraphErrors();
GraphErrors graph2q2 = new GraphErrors();
GraphErrors graph2q3 = new GraphErrors();


graph2q1.addPoint(12,carbon_differencebin4q1,0,deltacarbondiffmeanbin4q1);
graph2q2.addPoint(12,carbon_differencebin4q2,0,deltacarbondiffmeanbin4q2);
graph2q3.addPoint(12,carbon_differencebin4q3,0,deltacarbondiffmeanbin4q3);


graph2q1.addPoint(56,iron_differencebin4q1,0,deltairondiffmeanbin4q1);
graph2q2.addPoint(56,iron_differencebin4q2,0,deltairondiffmeanbin4q2);
graph2q3.addPoint(56,iron_differencebin4q3,0,deltairondiffmeanbin4q3);


graph2q1.addPoint(208,lead_differencebin4q1,0,deltaleaddiffmeanbin4q1);
graph2q2.addPoint(208,lead_differencebin4q2,0,deltaleaddiffmeanbin4q2);
graph2q3.addPoint(208,lead_differencebin4q3,0,deltaleaddiffmeanbin4q3);


graph2q1.setTitle("#Delta pT^2 vs Atomic Number for 4th interval of zh");
graph2q1.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph2q1.setTitleX("Atomic Number");


GraphErrors graph3q1 = new GraphErrors();
GraphErrors graph3q2 = new GraphErrors();
GraphErrors graph3q3 = new GraphErrors();


graph3q1.addPoint(12,carbon_differencebin1q1,0,deltacarbondiffmeanbin1q1);
graph3q2.addPoint(12,carbon_differencebin1q2,0,deltacarbondiffmeanbin1q2);
graph3q3.addPoint(12,carbon_differencebin1q3,0,deltacarbondiffmeanbin1q3);


graph3q1.addPoint(56,iron_differencebin1q1,0,deltairondiffmeanbin1q1);
graph3q2.addPoint(56,iron_differencebin1q2,0,deltairondiffmeanbin1q2);
graph3q3.addPoint(56,iron_differencebin1q3,0,deltairondiffmeanbin1q3);


graph3q1.addPoint(208,lead_differencebin1q1,0,deltaleaddiffmeanbin1q1);
graph3q2.addPoint(208,lead_differencebin1q2,0,deltaleaddiffmeanbin1q2);
graph3q3.addPoint(208,lead_differencebin1q3,0,deltaleaddiffmeanbin1q3);


graph3q1.setTitle("#Delta pT^2 vs Atomic Number for 1st interval of zh");
graph3q1.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph3q1.setTitleX("Atomic Number");


GraphErrors graph4q1 = new GraphErrors();
GraphErrors graph4q2 = new GraphErrors();
GraphErrors graph4q3 = new GraphErrors();


graph4q1.addPoint(12,carbon_differencebin5q1,0,deltacarbondiffmeanbin5q1);
graph4q2.addPoint(12,carbon_differencebin5q2,0,deltacarbondiffmeanbin5q2);
graph4q3.addPoint(12,carbon_differencebin5q3,0,deltacarbondiffmeanbin5q3);


graph4q1.addPoint(56,iron_differencebin5q1,0,deltairondiffmeanbin5q1);
graph4q2.addPoint(56,iron_differencebin5q2,0,deltairondiffmeanbin5q2);
graph4q3.addPoint(56,iron_differencebin5q3,0,deltairondiffmeanbin5q3);


graph4q1.addPoint(208,lead_differencebin5q1,0,deltaleaddiffmeanbin5q1);
graph4q2.addPoint(208,lead_differencebin5q2,0,deltaleaddiffmeanbin5q2);
graph4q3.addPoint(208,lead_differencebin5q3,0,deltaleaddiffmeanbin5q3);


graph4q1.setTitle("#Delta pT^2 vs Atomic Number for 5th interval of zh");
graph4q1.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph4q1.setTitleX("Atomic Number");


GraphErrors graph5q1 = new GraphErrors();
GraphErrors graph5q2 = new GraphErrors();
GraphErrors graph5q3 = new GraphErrors();


graph5q1.addPoint(12,carbon_differencebin6q1,0,deltacarbondiffmeanbin6q1);
graph5q2.addPoint(12,carbon_differencebin6q2,0,deltacarbondiffmeanbin6q2);
graph5q3.addPoint(12,carbon_differencebin6q3,0,deltacarbondiffmeanbin6q3);


graph5q1.addPoint(56,iron_differencebin6q1,0,deltairondiffmeanbin6q1);
graph5q2.addPoint(56,iron_differencebin6q2,0,deltairondiffmeanbin6q2);
graph5q3.addPoint(56,iron_differencebin6q3,0,deltairondiffmeanbin6q3);


graph5q1.addPoint(208,lead_differencebin6q1,0,deltaleaddiffmeanbin6q1);
graph5q2.addPoint(208,lead_differencebin6q2,0,deltaleaddiffmeanbin6q2);
graph5q3.addPoint(208,lead_differencebin6q3,0,deltaleaddiffmeanbin6q3);


graph5q1.setTitle("#Delta pT^2 vs Atomic Number for 6th interval of zh");
graph5q1.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph5q1.setTitleX("Atomic Number");



// I'm still not sure what this part with the directory does

String dirname = "/electron";
TDirectory dir = new TDirectory();
dir.mkdir(dirname);
dir.cd(dirname);


/*
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

*/


int c1a_title_size = 30;
TCanvas c1 = new TCanvas("q2",800,800);
// number columns number rows
c1.divide(3,2); 


c1.cd(0);
graph3q1.setMarkerColor(1);
graph3q1.setMarkerStyle(1);
graph3q2.setMarkerColor(2);
graph3q2.setMarkerStyle(2);
graph3q3.setMarkerColor(3);
graph3q3.setMarkerStyle(3);

c1.draw(graph3q1, "same");
c1.draw(graph3q2, "same");
c1.draw(graph3q3, "same");



// create 4 graphs of all of the q1, q2, etc.
// like: c1.draw(graph3q1, "same");
// .setMarkerColor()
// grpah.setMarkerSize()
// grpah.setMarkerStyle(1), etc.
// put this at the beginning right after the cd(1) stuff

c1.cd(1);
graphq1.setMarkerColor(1);
graphq1.setMarkerStyle(1);
graphq2.setMarkerColor(2);
graphq2.setMarkerStyle(2);
graphq3.setMarkerColor(3);
graphq3.setMarkerStyle(3);

c1.draw(graphq1, "same");
c1.draw(graphq2, "same");
c1.draw(graphq3, "same");



c1.cd(2);
graph1q1.setMarkerColor(1);
graph1q1.setMarkerStyle(1);
graph1q2.setMarkerColor(2);
graph1q2.setMarkerStyle(2);
graph1q3.setMarkerColor(3);
graph1q3.setMarkerStyle(3);

c1.draw(graph1q1, "same");
c1.draw(graph1q2, "same");
c1.draw(graph1q3, "same");




c1.cd(3);
graph2q1.setMarkerColor(1);
graph2q1.setMarkerStyle(1);
graph2q2.setMarkerColor(2);
graph2q2.setMarkerStyle(2);
graph2q3.setMarkerColor(3);
graph2q3.setMarkerStyle(3);

c1.draw(graph2q1, "same");
c1.draw(graph2q2, "same");
c1.draw(graph2q3, "same");




c1.cd(4);
graph4q1.setMarkerColor(1);
graph4q1.setMarkerStyle(1);
graph4q2.setMarkerColor(2);
graph4q2.setMarkerStyle(2);
graph4q3.setMarkerColor(3);
graph4q3.setMarkerStyle(3);

c1.draw(graph4q1, "same");
c1.draw(graph4q2, "same");
c1.draw(graph4q3, "same");




c1.cd(5);
graph5q1.setMarkerColor(1);
graph5q1.setMarkerStyle(1);
graph5q2.setMarkerColor(2);
graph5q2.setMarkerStyle(2);
graph5q3.setMarkerColor(3);
graph5q3.setMarkerStyle(3);

c1.draw(graph5q1, "same");
c1.draw(graph5q2, "same");
c1.draw(graph5q3, "same");



System.out.println("Black Square is the lowest range of q2");
System.out.println("Red Triangle is the middle range of q2");
System.out.println("Green Inverted Triangle is the highest range of q2");


