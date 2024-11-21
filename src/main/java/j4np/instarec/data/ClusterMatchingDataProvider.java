/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.data;

import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.instarec.utils.DataEntry;
import j4np.instarec.utils.DataList;
import j4np.utils.io.OptionParser;

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
public class ClusterMatchingDataProvider {
    
    public ClusterMatchingDataProvider(){
        
    }

    public Boolean isInAllLayers(Leaf part,int row){
      for(int i=18;i<27;i++){
        if(part.getDouble(i,row)==0){return false;}
        if(part.getDouble(29,row)==0){return false;}
      }
      return true;
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
    
    public void process(String file, String output, int limTotNegEvs){

      //used to convert ls to strips
      DataList LtoSconv = readCSV("LtoStrip_convTable.csv",69);

      //output training sample size per sector
      int[] count = new int[14];
      count[0]=0;
      count[7]=0;

      //remove existing output files so we can overwrite output
      delete_existing_output(output+"_negatives.csv");
      delete_existing_output(output+"_positives.csv");
      for(int i=1;i<7;i++){
        delete_existing_output(output+"_negatives_sector"+String.valueOf(i)+".csv");
        delete_existing_output(output+"_positives_sector"+String.valueOf(i)+".csv");
        count[i]=0;
        count[i+7]=0;
      }

      HipoReader r = new HipoReader(file);
      Event ev = new Event();
  
      Leaf part = new Leaf(42, 15, "i", 1200);
      
      while(r.hasNext() && count[0]<limTotNegEvs){
        r.nextEvent(ev);
        ev.read(part,42,15);
        for(int row=0;row<part.getRows();row++){
          //only take particles that hit FTOF and all layers of ECAL
          if(isInAllLayers(part,row)){

            StringBuilder csvLineBuilder = new StringBuilder();
            //add wires and ECAL Ls
            /*for (int j=12;j<27;j++){ 
              csvLineBuilder.append(String.format("%.6f,",part.getDouble(j,row)));
            }*/

            //best to conv to strips first
            //add wires
            for (int j=12;j<18;j++){ 
              csvLineBuilder.append(String.format("%.6f,",part.getDouble(j,row)));
            }
            
            float[] ls = new float[9];
            for (int j=18;j<27;j++){ 
              ls[j-18]=(float)part.getDouble(j,row);
            }
            
            int[] strips = new int[] { -2, -2, -2, -2, -2, -2, -2, -2, -2 };
            convLtoStrip(ls,strips,LtoSconv);

            for (int j=0;j<9;j++){ 
              csvLineBuilder.append(String.format("%d,",strips[j]));
            }

            csvLineBuilder.append(String.format("%.6f,",part.getDouble(29,row)));
            csvLineBuilder.append(String.format("%.6f",part.getDouble(27,row)));
            short charge = part.getShort(4,row);

            int countOffset=0;
            String chargeString="_negatives";
            if(charge==1){
              countOffset=7;
              chargeString="_positives";
            }

            //append string per event for all sectors
            try{
              count[0+countOffset]++;
              FileWriter writer = new FileWriter(output+chargeString+".csv",true);
              PrintWriter pw = new PrintWriter(writer);
              pw.println(csvLineBuilder.toString());
              pw.close();
            } catch (IOException e){
              e.printStackTrace();
            }

            try{
              count[part.getShort(3,row)+countOffset]++;
              FileWriter writer = new FileWriter(output+chargeString+"_sector"+String.valueOf(part.getShort(3,row))+".csv",true);
              PrintWriter pw = new PrintWriter(writer);
              pw.println(csvLineBuilder.toString());
              pw.close();
            } catch (IOException e){
              e.printStackTrace();
            }
    
          }
        }
      }

      System.out.println("\nRead file "+file);
      System.out.println("Training Data Statistics:");
      System.out.printf("\t Negatives: %d | Positives: %d \n",count[0],count[7]);
      for(int i=1;i<7;i++){
        System.out.printf("Sector %d Negatives: %d | Positives: %d \n",i,count[i],count[i+7]);
      }

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

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar wTrainCF.h5 -o training_data/cfTrain
    
    public static void main(String[] args){
        
      System.out.println("----- starting cluster finder data provider ");
      OptionParser p = new OptionParser();
      p.addRequired("-o", "output file name");
      p.parse(args);
        
      ClusterMatchingDataProvider dp = new ClusterMatchingDataProvider();
        
      dp.process(p.getInputList().get(0), p.getOption("-o").stringValue(),150000);
        
    }
}
