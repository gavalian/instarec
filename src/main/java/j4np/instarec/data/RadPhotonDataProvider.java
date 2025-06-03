/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.data;

import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.utils.io.OptionParser;
import twig.data.GraphErrors;
import twig.data.H1F;
import twig.data.H2F;
import twig.data.StatNumber;
import twig.data.TDirectory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

//for priting
import j4np.hipo5.data.Bank;

/**
 *
 * @author tyson
 */
public class RadPhotonDataProvider {
    
    public RadPhotonDataProvider(){
    }

    public void getPart(Bank RECPart, int pindex, double[] part){
      double pz = RECPart.getFloat("pz", pindex);
      double px = RECPart.getFloat("px", pindex);
      double py = RECPart.getFloat("py", pindex);
      double p=Math.sqrt(px*px+py*py+pz*pz);
      double Theta = Math.acos(pz / p)*(180/Math.PI);// Math.atan2(Math.sqrt(px*px+py*py),pz);
      double Phi = Math.atan2(py, px)*(180/Math.PI);
      int pid=RECPart.getInt("pid", pindex);
      part[0]=pid;
      part[1]=p;
      part[2]=Theta;
      part[3]=Phi;
    }

    public void getCalInfo(Bank RECCal, int pindex, double P, double Es[], double Ls[]){
    
      float PCALE=0,ECINE=0,ECOUTE=0;
      for (int k = 0; k < RECCal.getRows(); k++) {
        short i = RECCal.getShort("pindex", k);
        float lu=RECCal.getFloat("lu",k);
        float lv=RECCal.getFloat("lv",k);
        float lw=RECCal.getFloat("lw",k);
        float energy=RECCal.getFloat("energy",k);
        byte layer=RECCal.getByte("layer",k);
        if (i == pindex ) {
          //Cal_index=index;
          if(layer==1){
            PCALE=energy;
            Ls[layer-1]=lu;
            Ls[layer-1+1]=lv;
            Ls[layer-1+2]=lw;
          } else if(layer==4){
            ECINE=energy;
            Ls[layer-1]=lu;
            Ls[layer-1+1]=lv;
            Ls[layer-1+2]=lw;
          } else if(layer==7){
            ECOUTE=energy;
            Ls[layer-1]=lu;
            Ls[layer-1+1]=lv;
            Ls[layer-1+2]=lw; 
          }
        }
      }
      Es[0]=PCALE;
      Es[1]=ECINE;
      Es[2]=ECOUTE;
      Es[3]=(PCALE + ECINE + ECOUTE)/P;
      Es[4]=PCALE/P;
      Es[5]=ECINE/P;
      Es[6]=ECOUTE/P;
    }

    public Boolean getRadiatedPhoton(Bank RECPart, double[] electron, double photon[]){

      for (short i = 0; i < RECPart.getRows(); i++) {
        getPart(RECPart,i,photon);
        if(Math.abs(electron[2]-photon[2])<0.5 && photon[0]==22){
          return true;
        }
      }

      return false;

    }

    public Boolean passFid(double[] Ls){
      if(Ls[0]>9 && Ls[1]>9 && Ls[2]>9){
        return true;
      } else{
        return false;
      }
    }

    //returns false when all predicted ECAL strips or FTOF comp are zero
    //ie pred missed detector
    //this often indicates a bad track
    public Boolean trackOutOfDet(Leaf part, int row) {
  
      if (part.getDouble(18,row)==0 && part.getDouble(19,row)==0 && part.getDouble(20,row)==0) {
        return true;
      }
  
      if (part.getDouble(21,row)==0 && part.getDouble(22,row)==0 && part.getDouble(23,row)==0) {
        return true;
      }
  
      if (part.getDouble(24,row)==0 && part.getDouble(25,row)==0 && part.getDouble(26,row)==0) {
        return true;
      }

      if (part.getDouble(29,row)==0 ) {
        return true;
      }

      return false;
    }

    //returns false when all predicted ECAL strips or FTOF comp are zero
    //ie pred missed detector
    //this often indicates a bad track
    public int countTracksInSector(Leaf part, int sector) {

      int nSect=0;
      for(int row=0;row<part.getRows();row++){
        short sect = part.getShort(3,row);
        if(sect==sector){
          nSect++;
        }
      }
      return nSect;
    }
    
    public int process(String file, int limTotNegEvs,double threshold, String endName,String output){

      //remove existing output files so we can overwrite output
      delete_existing_output(output+".csv");
      //if did this per sector
      //for(int i=1;i<14;i++){
      //  if(i<7){delete_existing_output(output+"_sector"+String.valueOf(i)+".csv");}
      //}

      //output training sample size per sector
      int count = 0, countRad=0, countGoodInsta=0, countBadInsta=0;

      HipoReader r = new HipoReader(file);
      Event ev = new Event();
  
      Leaf part = new Leaf(32, 99, "i", 1200);
      Leaf pred_part = new Leaf(32, 3, "i", 1200);
      Leaf pred_ECAL = new Leaf(32, 1, "i", 1200);
      Leaf pred_HTCC = new Leaf(32, 98, "i", 1200);

      Bank recpart = r.getBank("REC::Particle");
      Bank reccal = r.getBank("REC::Calorimeter");


      double[] electron = new double[4];
      double[] elEs=new double[7];
      double[] elLs=new double[9];
      double[] photon = new double[4];
      double[] phEs=new double[7];
      double[] phLs=new double[9];

      //don't count if set limTotNegEvs=-1
      if(limTotNegEvs==-1){count=-2;}

      while(r.hasNext() && count<limTotNegEvs){

        r.nextEvent(ev);
        ev.read(part,32,99);
        ev.read(pred_part,32,3);
        ev.read(pred_HTCC,32,98);
        ev.read(pred_ECAL,32,1);
        ev.read(recpart);
        ev.read(reccal);
        

        //same rows in part bank and pred_part bank
        for(int row=0;row<pred_part.getRows();row++){
          short charge = pred_part.getShort(4,row);
          short sector = pred_part.getShort(3,row);
          int pindex=part.getShort(0,row);
          float resp=(float)pred_part.getDouble(2,row);

          double totHTCC=pred_part.getDouble(30,row)+pred_part.getDouble(31,row)+pred_part.getDouble(32,row);

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

          csvLineBuilder.append("0,1");

          //require neg parts in FD matched to rec track
          //apply fiducial cuts 
          if(sector>0 && charge==-1 && pindex!=999 ){
            for(int rowphoton=0;rowphoton<recpart.getRows();rowphoton++){ 

              getPart(recpart, pindex, electron);
              int recpid=(int)electron[0];
          
              getPart(recpart, rowphoton, photon);
              int recpidph=(int)photon[0];

              double dTheta=photon[2]-electron[2];
              double dPhi=photon[3]-electron[3];

              

              //different part for electron and photon
              if(rowphoton!=pindex && recpid!=11 && recpidph==22 && electron[1]>0.01){

                
                getCalInfo(reccal, pindex, electron[1], elEs, elLs);
                getCalInfo(reccal, rowphoton, photon[1], phEs, phLs);
        
                //make sure part is good electron candidate
                //fids && passFid(elLs)
                if(elEs[3]>0.01  && totHTCC!=0 && Math.abs(dTheta)<2.5){// elEs[0]>0.06 elEs[5] > (0.15 - elEs[4]) && elEs[3]>0.15

                  
                  double difLU=elLs[0]-phLs[0];
                  double difLV=elLs[1]-phLs[1];
                  double difLW=elLs[2]-phLs[2];
                  double dist=Math.sqrt(difLU*difLU+difLV*difLV+difLW*difLW);

                
                  //region with less bg from non rad photons
                  if( dPhi>-30 && dPhi<30){ //dist>15 &&
                    
                    if(Math.abs(dTheta)<0.5 ){
                      try{
                        FileWriter writer = new FileWriter(output+".csv",true);
                        PrintWriter pw = new PrintWriter(writer);
                        pw.println(csvLineBuilder.toString());
                        pw.close();
                        countRad++;
                        //don't count if set limTotNegEvs=-1
                        if(limTotNegEvs!=-1){
                          count++;
                        }
                      } catch (IOException e){
                        e.printStackTrace();
                      }
                      if(resp>threshold){
                        countGoodInsta++;
                      } else {
                        countBadInsta++;
                      }

                    }  
                  }
                }

              }
            }
          }
        }
      }

      float eff=(float)countGoodInsta/countRad;

      System.out.printf("\n\nNumber of missIDed electrons with radiated photon %d , recovered %d, not recovered %d\n",countRad,countGoodInsta,countBadInsta);
      System.out.printf("Recovery Efficiency %f\n",eff);

      return countRad;

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

    public void copyOverNegatives(int numLines, String inputPath, String outputPath) {
        try (
            BufferedReader reader = new BufferedReader(new FileReader(inputPath+".csv"));
            BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath+".csv", true)) // append = true
        ) {
            String line;
            int count = 0;

            while ((line = reader.readLine()) != null && count < numLines) {
                if (line.trim().endsWith("1,0")) {
                    writer.write(line);
                    writer.newLine();
                    count++;
                }
            }
            System.out.println("Appended " + count + " negatives to: " + outputPath);
        } catch (IOException e) {
            System.err.println("Error processing files: " + e.getMessage());
        }
    }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar -in w.h5 -o training_data/ElPIDTrainRadPhotons -nr training_data/ElPIDTrain
    //read plots in j4shell with eg TwigStudio.browser("plots/RadPhoton0.twig");
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting radiated photon validator ");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.addRequired("-o", "output file name");
      p.addRequired("-nr", "usual (not rad) training file name");
      p.parse(args);

      String endName="_noFid"; // used to change output path of plots (eg adding _NoFiducialCuts)
        
      RadPhotonDataProvider dp = new RadPhotonDataProvider();
      int nPositives=dp.process(p.getOption("-in").stringValue(),-1,0.05,endName,p.getOption("-o").stringValue());
      dp.copyOverNegatives(nPositives,p.getOption("-nr").stringValue(),p.getOption("-o").stringValue());

  
      
      
        
    }
}
