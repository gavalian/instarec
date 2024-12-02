/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.networks;

import deepnetts.data.TabularDataSet;
import deepnetts.net.FeedForwardNetwork;
import j4np.instarec.utils.DataEntry;
import j4np.instarec.utils.DataList;
import j4np.instarec.utils.DeepNettsIO;
import j4np.instarec.utils.DeepNettsTrainer;
import j4np.instarec.utils.EntryTransformer;
import j4np.instarec.utils.NeuralModel;
import j4np.utils.io.OptionParser;
import j4np.utils.io.TextFileReader;
import twig.data.H1F;
import twig.data.TDirectory;
import twig.studio.TwigStudio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
import javax.visrec.ml.data.DataSet;


/**
 *
 * @author tyson
 */
public class TrainingClusterFinder {

    static int nPredictedOutput=11;

    //Will get rid of this, there seems to be a bug in DataList.csv() method
    public static DataList readCSV(String file, int nIn, int nOut){
      DataList list = new DataList();

      try {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        while ((line = reader.readLine()) != null) {

          String[] tokens = line.split(",");
          if (tokens.length != nIn + nOut) {
            throw new IllegalArgumentException("Invalid line length: " + line);
          }

          // Create arrays for nIn and nOut
          float[] inputArray = new float[nIn];
          float[] outputArray = new float[nOut];

          // Fill inputArray
          for (int i = 0; i < nIn; i++) {
            inputArray[i] = Float.parseFloat(tokens[i].trim());
          }

          // Fill outputArray
          for (int i = 0; i < nOut; i++) {
              outputArray[i] = Float.parseFloat(tokens[nIn + i].trim());
          }

         list.add(new DataEntry(inputArray, outputArray));

        }

        reader.close();

      } catch (IOException e) {
          System.err.println("Error reading file: " + e.getMessage());
      } catch (NumberFormatException e) {
          System.err.println("Error parsing number: " + e.getMessage());
      }

      return list;
    }
    
    
    public static void trainNetwork(String dataPath,String networkPath){
        
      FeedForwardNetwork network = DeepNettsTrainer.createRegression(new int[]{6, 12, 36, 27, 18,nPredictedOutput});
      EntryTransformer     trans = new EntryTransformer();
        
      trans.input().add(6,1, 112); //wires
      trans.output().add(9, 1, 500); //ECAL locs
      trans.output().add(1, 62); //FTOF component
      trans.output().add(1, 1000); //FTOF path
        
      //DataList data = DataList.csv(dataPath + ".csv",
      // DataList.range(0, 6), DataList.range(6, 6 + nPredictedOutput)); //_outbending

      DataList data=readCSV(dataPath+".csv",6,nPredictedOutput); 

      data.shuffle();
      data.scan();
      //data.show();
        
      DeepNettsTrainer   trainer = new DeepNettsTrainer(network);
      trainer.updateLayers();        
      DataSet dset = trainer.convert(data, trans);        
      trainer.train(dset, 1000);
      trainer.save(networkPath+".json", trans);
    }
    
    public static void testNetwork(String dataPath,String networkPath, int sector,String charge){
      String networkPath_ws=networkPath;
      String dataPath_ws=dataPath;
      if(sector!=0){
        networkPath_ws=networkPath+"_sector"+String.valueOf(sector);
        dataPath_ws=dataPath+"_sector"+String.valueOf(sector);
      }
      NeuralModel model = NeuralModel.jsonFile(networkPath_ws+".json");
      //System.out.println(model.summary());
        
      //DataList data = DataList.csv(dataPath_ws + ".csv",
      //  DataList.range(0, 6), DataList.range(6, 6 + nPredictedOutput)); //_outbending
      DataList data=readCSV(dataPath_ws+".csv",6,nPredictedOutput); 
      data.shuffle();

      String[] calLayer = new String[3];
      calLayer[0]="PCAL";
      calLayer[1]="ECIN";
      calLayer[2]="ECOUT";
        
      for(int i=0;i<3;i++){

        int out_start=3*i;
        float[] predicted = new float[nPredictedOutput];

        H1F hDifU = new H1F("U View ( "+calLayer[i]+", "+charge+", "+String.valueOf(sector)+")", 21, -10.5, 10.5);
        hDifU.attr().setLineColor(2);
        hDifU.attr().setLineWidth(3);
        hDifU.attr().setTitleX("Position Difference [strips]");

        H1F hDifV = new H1F("V View ( "+calLayer[i]+", "+charge+", "+String.valueOf(sector)+")", 21, -10.5, 10.5);
        hDifV.attr().setLineColor(5);
        hDifV.attr().setLineWidth(3);
        hDifV.attr().setTitleX("Position Difference [strips]");

        H1F hDifW = new H1F("W View ( "+calLayer[i]+", "+charge+", "+String.valueOf(sector)+")", 21, -10.5, 10.5);
        hDifW.attr().setLineColor(3);
        hDifW.attr().setLineWidth(3);
        hDifW.attr().setTitleX("Position Difference [strips]");

        H1F hDifP = new H1F("FTOF Path ("+charge+" "+String.valueOf(sector)+")", 100, -10, 10);
        hDifP.attr().setLineColor(2);
        hDifP.attr().setLineWidth(3);
        hDifP.attr().setTitleX("Path Difference [cm]");

        H1F hDif = new H1F("FTOF Component ("+charge+" "+String.valueOf(sector)+")", 11, -5.5, 5.5);
        hDif.attr().setLineColor(2);
        hDif.attr().setLineWidth(3);
        hDif.attr().setTitleX("Position Difference [component]");

        for(int j = 0; j < data.getList().size(); j++){
          model.predict(data.getList().get(j).features(), predicted);
          float[] desired = data.getList().get(j).floatSecond();
          if(j<10 && i==0){
            if(j<1){System.out.println("Testing -- first 10 entries");}
            System.out.println("PCAL U/V/W ECIN U/V/W ECOUT U/V/W FTOF Commp/Path");
            for(int k=0;k<10;k++){System.out.printf("%d, ",Math.round(desired[k]));}
            System.out.printf("%d \n",Math.round(desired[10]));
            for(int k=0;k<10;k++){System.out.printf("%d, ",Math.round(predicted[k]));}
            System.out.printf("%d \n\n",Math.round(predicted[10]));
          }
          
          hDifU.fill((desired[out_start + 0] - predicted[out_start + 0]));
          hDifV.fill((desired[out_start + 1] - predicted[out_start + 1]));
          hDifW.fill((desired[out_start + 2] - predicted[out_start + 2]));
          hDif.fill((desired[9] - predicted[9]));
          hDifP.fill((desired[10] - predicted[10]));
        }
        TDirectory.export("plots/clusterfinder"+String.valueOf(sector)+".twig","/ai/training/"+String.valueOf(calLayer[i])+"/"+charge,hDifU);
        TDirectory.export("plots/clusterfinder"+String.valueOf(sector)+".twig","/ai/training/"+String.valueOf(calLayer[i])+"/"+charge,hDifV);
        TDirectory.export("plots/clusterfinder"+String.valueOf(sector)+".twig","/ai/training/"+String.valueOf(calLayer[i])+"/"+charge,hDifW);
        if(i==0){
          TDirectory.export("plots/clusterfinder"+String.valueOf(sector)+".twig","/ai/training/FTOF/"+charge,hDif);
          TDirectory.export("plots/clusterfinder"+String.valueOf(sector)+".twig","/ai/training/FTOF/"+charge,hDifP);
        }
      }
    }
    
    public static void transferTraining(String dataPath,String networkPath, int sector){
      String networkPath_ws=networkPath;
      String dataPath_ws=dataPath;
      if(sector!=0){
        networkPath_ws=networkPath+"_sector"+String.valueOf(sector);
        dataPath_ws=dataPath+"_sector"+String.valueOf(sector);
      }
      FeedForwardNetwork network = DeepNettsIO.read(networkPath+".json");     
      EntryTransformer   transformer = DeepNettsIO.getTransformer(networkPath+".json");
        
      DeepNettsTrainer   trainer = new DeepNettsTrainer(network);
      trainer.updateLayers();
        
      //DataList data = DataList.csv(dataPath_ws + ".csv",
      //  DataList.range(0, 6), DataList.range(6, 6 + nPredictedOutput)); //_outbending
      DataList data=readCSV(dataPath_ws+".csv",6,nPredictedOutput); 
      data.shuffle();
      data.scan();
      //data.show();
        
      DataSet dset = trainer.convert(data, transformer);        
      trainer.train(dset, 250);
      trainer.save(networkPath_ws+".json", transformer);
    }
    
    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar -in training_data/cfTrain -name etc/networks/clusterfinder/cf
    //read plots in j4shell with eg TwigStudio.browser("plots/clusterfinder0.twig");

    public static void main(String[] args){

      System.out.println("----- starting cluster finder training");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.addRequired("-name", "network name");
      p.parse(args);
        
      String[] charge = new String[2];
      charge[0]="negatives";
      charge[1]="positives";

      for(int i=0;i<2;i++){
        System.out.println("\n\n\nTraining for "+charge[i]);
        String dataPath=p.getOption("-in").stringValue()+"_"+charge[i];
        String networkPath=p.getOption("-name").stringValue()+"_"+charge[i];
        TrainingClusterFinder.trainNetwork(dataPath,networkPath);
        TrainingClusterFinder.testNetwork(dataPath,networkPath,0,charge[i]);
        for(int j=1;j<7;j++){
          System.out.printf("\nTransfer training for sector %d\n",j);
          TrainingClusterFinder.transferTraining(dataPath,networkPath,j);
          TrainingClusterFinder.testNetwork(dataPath,networkPath,j,charge[i]);
        }
      }
      
    }
}
