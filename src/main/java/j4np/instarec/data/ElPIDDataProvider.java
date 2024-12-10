/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.data;

import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Bank;
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
public class ElPIDDataProvider {
    
    public ElPIDDataProvider(){
        
    }

    public Boolean radiatedPhoton(Bank RECPart, int pindex){

      double pz = RECPart.getFloat("pz", pindex);
      double px = RECPart.getFloat("px", pindex);
      double py = RECPart.getFloat("py", pindex);
      double p=Math.sqrt(px*px+py*py+pz*pz);
      double Theta = Math.acos(pz / p)*(180/Math.PI);// Math.atan2(Math.sqrt(px*px+py*py),pz);
      double Phi = Math.atan2(py, px)*(180/Math.PI);
      for (short i = 0; i < RECPart.getRows(); i++) {
        int charge = RECPart.getByte("charge", i);
        double pzphoton = RECPart.getFloat("pz", i);
        double pxphoton = RECPart.getFloat("px", i);
        double pyphoton = RECPart.getFloat("py", i);
        double pphoton=Math.sqrt(pxphoton*pxphoton+pyphoton*pyphoton+pzphoton*pzphoton);
        double Thetaphoton = Math.acos(pzphoton / p)*(180/Math.PI);// Math.atan2(Math.sqrt(pxphoton*pxphoton+pyphoton*pyphoton),pzphoton);
        double Phiphoton = Math.atan2(pyphoton, pxphoton)*(180/Math.PI);
        //just use charge and not pid because some photons missided as neutrons
        if(Math.abs(Theta-Thetaphoton)<0.5 && charge==0){
          return true;
        }
      }

      return false;

    }

    
    public void process(String file, String output, int limTotNegEvs){

      //output training sample size per sector
      int[] count = new int[14];
      count[0]=0;

      //remove existing output files so we can overwrite output
      delete_existing_output(output+".csv");
      for(int i=1;i<14;i++){
        if(i<7){delete_existing_output(output+"_sector"+String.valueOf(i)+".csv");}
        count[i]=0;
      }

      HipoReader r = new HipoReader(file);
      Event ev = new Event();
  
      Leaf part = new Leaf(32200, 99, "i", 1200);
      Leaf pred_ECAL = new Leaf(32200, 2, "i", 1200);
      Leaf pred_HTCC = new Leaf(32200, 98, "i", 1200);

      Bank recpart = r.getBank("REC::Particle");

      int badEl=0;
      
      while(r.hasNext() && count[0]<limTotNegEvs){
        r.nextEvent(ev);
        ev.read(part,32200, 99);
        ev.read(pred_ECAL,32200,2);
        ev.read(pred_HTCC,32200,98);
        ev.read(recpart);
        for(int row=0;row<part.getRows();row++){

          short charge = part.getShort(4,row);
          short sector = part.getShort(3,row);
          int pid=part.getInt(1,row);
          int pindex=part.getShort(0,row);
          //only take negative particles 
          //and tracks matched to rec particles
          if(charge==-1 && sector>0 && pindex!=999){

            StringBuilder csvLineBuilder = new StringBuilder(); 

            float[] pred_strips = new float[9];
            float[] energy= new float[3];
            float sumE=0;
            for(int j=0;j<3;j++){energy[j]=0;}
            float[] DUs = new float[9];
            for(int j=0;j<pred_ECAL.getRows();j++){
              //match cluster to track
              if(pred_ECAL.getShort(0,j)==row){
                if(pred_ECAL.getInt(2,j)!=0){
                  pred_strips[pred_ECAL.getInt(2,j)-1]=(float)pred_ECAL.getDouble(3,j);
                  DUs[pred_ECAL.getInt(2,j)-1]=pred_ECAL.getInt(5,j);
                  sumE+=pred_ECAL.getDouble(4,j);
                  //sum energy in PCAL, ECIN, ECOUT
                  if(pred_ECAL.getInt(2,j)<4){
                    energy[0]+=pred_ECAL.getDouble(4,j);
                  } else if(pred_ECAL.getInt(2,j)>=4 && pred_ECAL.getInt(2,j)<7){
                    energy[1]+=pred_ECAL.getDouble(4,j);
                  } else if(pred_ECAL.getInt(2,j)>=7){
                    energy[2]+=pred_ECAL.getDouble(4,j);
                  }
                }
              }
            }

            csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,",energy[0],energy[1],energy[2]));

            for(int j=0;j<9;j++){
              csvLineBuilder.append(String.format("%.6f,",pred_strips[j]));
            }
            for(int j=0;j<9;j++){
              csvLineBuilder.append(String.format("%.6f,",DUs[j]));
            }

            //add wires
            for (int j=12;j<18;j++){ 
              csvLineBuilder.append(String.format("%.6f,",part.getDouble(j,row)));
            }
            //add HTCC, has one row per sector
            float sumHTCC=0;
            for(int j=0;j<8;j++){
              csvLineBuilder.append(String.format("%.6f,",pred_HTCC.getDouble(j+1,sector-1)));
              sumHTCC+=pred_HTCC.getDouble(j+1,sector-1);
            }

            int sectorBefore=sector-1;
            int sectorAfter=sector+1;
            if(sector==1){sectorBefore=6;}
            else if(sector==6){sectorAfter=1;}

            for(int j=0;j<8;j++){
              csvLineBuilder.append(String.format("%.6f,",pred_HTCC.getDouble(j+1,sectorBefore-1)));
              sumHTCC+=pred_HTCC.getDouble(j+1,sectorBefore-1);
            }

            for(int j=0;j<8;j++){
              csvLineBuilder.append(String.format("%.6f,",pred_HTCC.getDouble(j+1,sectorAfter-1)));
              sumHTCC+=pred_HTCC.getDouble(j+1,sectorAfter-1);
            }

            int countOffset=0,write=0,writesect=0;
            if(pid==11){
              csvLineBuilder.append("0,1");
              if(sumHTCC==0){badEl++;}
              //require at least some energy in cal and hit in HTCC to avoid bad tracks
              if(sumE!=0 && sumHTCC!=0){
                if(count[0]<count[7]){
                  write=1;
                }
                if(count[sector]<count[7+sector] ){
                  writesect=1;
                }
                //write=1;
                //writesect=1;
              }
            } else {
              countOffset=7;
              csvLineBuilder.append("1,0");
              //require at least some energy in cal to avoid bad tracks
              //radiated photons are a clear indication that the particle is an e-
              if(sumE!=0 && !radiatedPhoton(recpart, pindex)){
                /*if(count[7]<count[0]){
                  write=1;
                }
                if(count[7+sector]<count[sector] ){
                  writesect=1;
                }*/
                write=1;
                writesect=1;
              }
            }

            //balance dataset by only writing bg if have less than signal
            //assumes more bg than signal in data (usually true)
            //(if not true, swap equalities and write stuff in block above)
            if(write==1){
              try{
                count[0+countOffset]++;
                FileWriter writer = new FileWriter(output+".csv",true);
                PrintWriter pw = new PrintWriter(writer);
                pw.println(csvLineBuilder.toString());
                pw.close();
              } catch (IOException e){
                e.printStackTrace();
              }
            }

            if(writesect==1){
              try{
                count[sector+countOffset]++;
                FileWriter writer = new FileWriter(output+"_sector"+String.valueOf(sector)+".csv",true);
                PrintWriter pw = new PrintWriter(writer);
                pw.println(csvLineBuilder.toString());
                pw.close();
              } catch (IOException e){
                e.printStackTrace();
              }
            }
    
          }
        }
      }

      System.out.println("\nRead file "+file);
      System.out.println("Training Data Statistics:");
      System.out.printf("All sectors: e- %d other %d \n",count[0],count[7]);
      for(int i=1;i<7;i++){
        System.out.printf("Sector %d : e- %d other %d \n",i,count[i],count[i+7]);
      }
      System.out.printf("e- with no HTCC hits %d \n",badEl);

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

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar wTrainElPID.h5 -o training_data/ElPIDTrain
    
    public static void main(String[] args){
        
      System.out.println("----- starting electron PID data provider ");
      OptionParser p = new OptionParser();
      p.addRequired("-o", "output file name");
      p.parse(args);
        
      ElPIDDataProvider dp = new ElPIDDataProvider();
        
      dp.process(p.getInputList().get(0), p.getOption("-o").stringValue(),150000);
        
    }
}
