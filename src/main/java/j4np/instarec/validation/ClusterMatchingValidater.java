/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.validation;

import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.instarec.networks.TrainingClusterFinder;
import j4np.instarec.utils.DataEntry;
import j4np.instarec.utils.DataList;
import j4np.instarec.utils.NeuralModel;
import j4np.utils.io.OptionParser;
import twig.data.H1F;
import twig.data.TDirectory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 *
 * @author tyson
 */
public class ClusterMatchingValidater {
    
    public ClusterMatchingValidater(){
    }

    // convert local U/V/W position in ECAL to strip numbers
    // need strip numbers to read ECAL ADC bank
    void convLtoStrip(float[] Ls, int[] strips, DataList LtoSconv) {

      int[] nStrips = new int[] { 68, 63, 63, 37, 37, 37, 37, 37, 37 };
      for (int i = 0; i < 9; i++) {
        Ls[i] = Ls[i];
        float[] conv_det_view = LtoSconv.getList().get(i).floatFirst();
        for (int str = 0; str < (nStrips[i]); str++) {
          // in PCAL V, W take distance from other end
          // System.out.printf("wtf %d %d \n",i,str);
          if (i == 1 || i == 2) {
            if (Ls[i] < conv_det_view[str] && Ls[i] >= conv_det_view[str + 1]) {
              strips[i] = str + 1; // add 1 as the upper limit applies to strip+1 due to inverted order
            }
          } else {
            if (Ls[i] >= conv_det_view[str] && Ls[i] < conv_det_view[str + 1]) {
              strips[i] = str;
            }
          }
        }
      }
    }

    //Will get rid of this, there seems to be a bug in DataList.csv() method
    public static DataList readCSV(String file, int nIn){
      DataList list = new DataList();

      try {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        int lineNb=0;
        while ((line = reader.readLine()) != null) {

          String[] tokens = line.split(",");
          if (tokens.length != nIn ) {
            throw new IllegalArgumentException("Invalid line length: " + line);
          }

          // Create arrays for nIn and nOut
          float[] inputArray = new float[nIn];

          // Fill inputArray
          for (int i = 0; i < nIn; i++) {
            inputArray[i] = Float.parseFloat(tokens[i].trim());
          }

          float[] outputArray = new float[1];
          outputArray[0]=lineNb;
          list.add(new DataEntry(inputArray, outputArray));
          lineNb++;

        }

        reader.close();

      } catch (IOException e) {
          System.err.println("Error reading file: " + e.getMessage());
      } catch (NumberFormatException e) {
          System.err.println("Error parsing number: " + e.getMessage());
      }

      return list;
    }

    
    public void process(String file, int limTotNegEvs, short desired_sector, String desired_chargeSt, int calLayer){

      //used to convert ls to strips
      DataList LtoSconv = readCSV("LtoStrip_convTable.csv",69);

      //output training sample size per sector
      int count = 0;

      HipoReader r = new HipoReader(file);
      Event ev = new Event();
  
      Leaf part = new Leaf(42, 15, "i", 1200);
      Leaf pred_ECAL = new Leaf(42, 12, "i", 1200);
      Leaf pred_FTOF = new Leaf(42, 14, "i", 1200);

      
      String calLayerSt="PCAL";
      if(calLayer==1){
        calLayerSt="ECIN";
      } else if(calLayer==2){
        calLayerSt="ECOUT";
      }

      int desired_charge=-1;
      if(desired_chargeSt=="positives"){desired_charge=1;}

      H1F hDifU = new H1F("U View ( "+calLayerSt+", "+desired_chargeSt+", "+String.valueOf(desired_sector)+")", 21, -10.5, 10.5);
      hDifU.attr().setLineColor(2);
      hDifU.attr().setLineWidth(3);
      hDifU.attr().setTitleX("Position Difference [strips]");

      H1F hDifV = new H1F("V View ( "+calLayerSt+", "+desired_chargeSt+", "+String.valueOf(desired_sector)+")", 21, -10.5, 10.5);
      hDifV.attr().setLineColor(5);
      hDifV.attr().setLineWidth(3);
      hDifV.attr().setTitleX("Position Difference [strips]");

      H1F hDifW = new H1F("W View ( "+calLayerSt+", "+desired_chargeSt+", "+String.valueOf(desired_sector)+")", 21, -10.5, 10.5);
      hDifW.attr().setLineColor(3);
      hDifW.attr().setLineWidth(3);
      hDifW.attr().setTitleX("Position Difference [strips]");

      H1F hDifP = new H1F("FTOF Path ("+desired_chargeSt+" "+String.valueOf(desired_sector)+")", 100, -10, 10);
      hDifP.attr().setLineColor(2);
      hDifP.attr().setLineWidth(3);
      hDifP.attr().setTitleX("Path Difference [cm]");

      H1F hDif = new H1F("FTOF Component ("+desired_chargeSt+" "+String.valueOf(desired_sector)+")", 11, -5.5, 5.5);
      hDif.attr().setLineColor(2);
      hDif.attr().setLineWidth(3);
      hDif.attr().setTitleX("Position Difference [component]");
      
      while(r.hasNext() && count<limTotNegEvs){

        int out_start=3*calLayer;

        r.nextEvent(ev);
        ev.read(part,42,15);
        ev.read(pred_ECAL,42,12);
        ev.read(pred_FTOF,42,14);
        for(int row=0;row<part.getRows();row++){


          float[] ls = new float[9];
          for (int j=18;j<27;j++){ 
            ls[j-18]=(float)part.getDouble(j,row);
          }

          float[] track = new float[6];
          for (int j=0;j<6;j++){ 
            track[j]=(float)part.getDouble(j+12,row);
          }
          
          int[] strips = new int[] { -2, -2, -2, -2, -2, -2, -2, -2, -2,-2,-2 };
          convLtoStrip(ls,strips,LtoSconv);
          strips[9]=(int)Math.round(part.getDouble(29,row));
          strips[10]=(int)Math.round(part.getDouble(27,row));

          float[] pred_strips = new float[11];
          for(int j=0;j<pred_ECAL.getRows();j++){
            //match cluster to track
            if(pred_ECAL.getShort(0,j)==row){
              if(pred_ECAL.getInt(2,j)!=0){
                pred_strips[pred_ECAL.getInt(2,j)-1]=(float)pred_ECAL.getDouble(3,j);
              }
            }
          }
          for(int j=0;j<pred_FTOF.getRows();j++){
            //match cluster to track, take layer 2
            if(pred_FTOF.getShort(0,j)==row && pred_FTOF.getInt(2,j)==2){
              pred_strips[9]=(float)pred_FTOF.getDouble(3,j); //pred comp
              pred_strips[10]=(float)pred_FTOF.getDouble(4,j); //pred path
            }
          }

          short charge = part.getShort(4,row);
          short sector=part.getShort(3,row);

          Boolean strip_nn=true;
          Boolean pred_strip_nn=true;
          for(int j=out_start;j<out_start+3;j++){
            if(strips[j]<1){strip_nn=false;}
            if(pred_strips[j]<1){pred_strip_nn=false;}
            if(out_start==0){
              //for pcal layer also plot FTOF
              if(strips[9]<1 || strips[10]<1){strip_nn=false;}
              if(pred_strips[9]<1 || pred_strips[10]<1){pred_strip_nn=false;}
            }
            //System.out.printf("Lay %d , true %d pred %f\n",j,strips[j],pred_strips[j]);
          }

          //plot for all when desired_sector is 0
          if(desired_sector==0){sector=desired_sector;}

          //only plot when a layer is not all 0
          //should probably also check what happens when layers are 0
          //does the prediction find non null sector?
          if(charge==desired_charge && sector==desired_sector && strip_nn && pred_strip_nn){
            hDifU.fill((strips[out_start + 0] - pred_strips[out_start + 0])); 
            hDifV.fill((strips[out_start + 1] - pred_strips[out_start + 1]));
            hDifW.fill((strips[out_start + 2] - pred_strips[out_start + 2]));
            hDif.fill((strips[9] - pred_strips[9]));
            hDifP.fill((strips[10] - pred_strips[10]));
            count++;
          }
          
        }
      }
      TDirectory.export("plots/clusterfinder"+String.valueOf(desired_sector)+".twig","/ai/validation/"+String.valueOf(calLayerSt)+"/"+desired_chargeSt,hDifU);
      TDirectory.export("plots/clusterfinder"+String.valueOf(desired_sector)+".twig","/ai/validation/"+String.valueOf(calLayerSt)+"/"+desired_chargeSt,hDifV);
      TDirectory.export("plots/clusterfinder"+String.valueOf(desired_sector)+".twig","/ai/validation/"+String.valueOf(calLayerSt)+"/"+desired_chargeSt,hDifW);
      if(calLayer==0){
        TDirectory.export("plots/clusterfinder"+String.valueOf(desired_sector)+".twig","/ai/validation/FTOF/"+desired_chargeSt,hDif);
        TDirectory.export("plots/clusterfinder"+String.valueOf(desired_sector)+".twig","/ai/validation/FTOF/"+desired_chargeSt,hDifP);
      }

      System.out.println("\n\nValidator Statistics:");
      System.out.printf(calLayerSt+" Sector %d, Charge %d\n",desired_sector,desired_charge);
      System.out.printf("N Events: %d \n",count);
      
    }

    
    public void delete_existing_output(String output){
      // Specify the path where you want to save the CSV file
      Path filePath = Paths.get(output);

      // Delete the file if it already exists
      try {
        Files.deleteIfExists(filePath);
      } catch (IOException e) {
        e.printStackTrace();
      }
    }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar -in w.h5 
    //read plots in j4shell with eg TwigStudio.browser("plots/clusterfinder0.twig");
    
    public static void main(String[] args){
        
      System.out.println("----- starting cluster finder data validater ");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.parse(args);
        
      ClusterMatchingValidater dp = new ClusterMatchingValidater();

      String[] chargeSt = new String[2];
      chargeSt[0]="negatives";
      chargeSt[1]="positives";

      //repeating reading of file 3*2*7 times
      // ==> innefficient
      //doing so to avoid writing code for 3*3*2*7 + 2*2*7 histos ...
      for(int lay=0;lay<3;lay++){
        for(int charge=0;charge<2;charge++){
          for(short sector=0;sector<1;sector++){ //7
            dp.process(p.getOption("-in").stringValue(),150000,sector,chargeSt[charge],lay);
          }
        }
      }
      
        
    }
}
