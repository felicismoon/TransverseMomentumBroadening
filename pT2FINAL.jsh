// this creates histograms of the pT^2 data for the deuterium and heavy targets for the 3 elements. It also creates a graph of delta pT^2 vs. A

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


// note: the deuterium target is iTgt==0, while the heavy target is iTgt==1
List vec = tree1.getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt==0",entries1);
List vec1 = tree1.getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt==1",entries1);
List vec3 = tree2.getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt==1",entries2);
List vec2 = tree2.getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt==0",entries2); 
List vec5 = tree3.getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt==1",entries3);
List vec4 = tree3.getDataVectors("pT2","pFidCut==1&&eFidCut==1&&iTgt==0",entries3); 




H1F[] histyArray = new H1F[6];


String[] hTitle = {"pT^2 vs Counts (Deuterium and Carbon)", "pT^2 vs Counts (Heavy Target and Carbon)", "pT^2 vs Counts (Deuterium and Iron)", "pT^2 vs Counts (Heavy Target and Iron)", "pT^2 vs Counts (Deuterium and Lead)", "pT^2 vs Counts (Heavy Target and Lead)"};

List[] vecArray = new List[6];
vecArray[0] = vec;
vecArray[1] = vec1;
vecArray[2] = vec2;
vecArray[3] = vec3;
vecArray[4] = vec4;
vecArray[5] = vec5;



for (int i = 0; i < 6; i++) {	
	String hname = "hpT2_" + i;
	histyArray[i] = new H1F().create(hname, 100, (DataVector)vecArray[i].get(0),0,1);
	histyArray[i].setTitle(hTitle[i]);
}



for (int i = 0; i < 6; i++) {	
	histyArray[i].setFillColor(i+2);
	histyArray[i].setTitleX("pT^2 (GeV^2/c^2)");
	histyArray[i].setTitleY("Counts");
}


// this part determines the average value of pT^2 for each histogram
double[] histyMean = new double[6];
for (int i = 0; i < 6; i++) {   
	histyMean[i] = histyArray[i].getMean();
}


double mean, mean1, mean2, mean3, mean4, mean5;
mean = histyMean[0];
mean1 = histyMean[1];
mean2 = histyMean[2];
mean3 = histyMean[3];
mean4 = histyMean[4];
mean5 = histyMean[5];


// this part determines the values of delta pT^2 for each element
double carbon_difference = mean1-mean; 
double iron_difference = mean2-mean3; 
double lead_difference = mean4-mean5; 


// this part calculates the error bars for the graph
double[] errorMean = new double[6];
for (int i = 0; i < 6; i++) {	
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



double deltacarbondiffmean = 0;
double deltairondiffmean = 0;
double deltaleaddiffmean = 0;


deltacarbondiffmean = Math.sqrt((errorMean[0])*(errorMean[0])-(errorMean[1])*(errorMean[1]));
deltairondiffmean = Math.sqrt((errorMean[2])*(errorMean[2])-(errorMean[3])*(errorMean[3]));
deltaleaddiffmean = Math.sqrt((errorMean[4])*(errorMean[4])-(errorMean[5])*(errorMean[5])); 



GraphErrors graph = new GraphErrors();
graph.addPoint(12,carbon_difference,0,deltacarbondiffmean);
graph.addPoint(56,iron_difference,0,deltairondiffmean);
graph.addPoint(208,lead_difference,0,deltaleaddiffmean);
graph.setTitle("#Delta pT^2 vs Atomic Number");
graph.setTitleY("#Delta pT^2 (GeV^2/c^2)");
graph.setTitleX("Atomic Number");


// I'm still not sure what this part with the directory does
String dirname = "/electron";
TDirectory dir = new TDirectory();
dir.mkdir(dirname);
dir.cd(dirname);

F1D func = new F1D("func","[a]+[b]*x", 0, 210);
DataFitter.fit(func, graph, "Q");
F1D func1 = new F1D("func1","[a]+[b]*x+[c]*x*x", 0, 210);
DataFitter.fit(func1, graph, "Q");
func.setLineColor(6);
func.setLineWidth(3);
func1.setLineColor(7);
func1.setLineWidth(3);



int c1a_title_size = 24;
TCanvas c1 = new TCanvas("c1",800,800);
c1.divide(2,4); 


for (int i = 0; i < 6; i++) {
	c1.cd(i);
	c1.draw(histyArray[i]);
	dir.addDataSet(histyArray[i]);
}



c1.cd(6);
c1.draw(graph);
c1.draw(func, "same");
c1.draw(func1, "same");




// this code creates 6 hists if you want to comment out the pt graph

/*c1.cd(0);
c1.draw(histyArray[0]);

c1.cd(1);
c1.draw(histyArray[2]);

c1.cd(2);
c1.draw(histyArray[4]);

c1.cd(3);
c1.draw(histyArray[1]);

c1.cd(4);
c1.draw(histyArray[3]);

c1.cd(5);
c1.draw(histyArray[5]); 
*/




System.out.println("Deuterium and Carbon Average: "+mean);
System.out.println("Heavy and Carbon Average: "+mean1);
System.out.println("Deuterium and Iron Average: "+mean2);
System.out.println("Heavy and Iron Average: "+mean3);
System.out.println("Deuterium and Lead Average: "+mean4);
System.out.println("Heavy and Lead Average: "+mean5);
System.out.println("Carbon Difference: "+carbon_difference);
System.out.println("Iron Difference: "+iron_difference);
System.out.println("Lead Difference: "+lead_difference);





