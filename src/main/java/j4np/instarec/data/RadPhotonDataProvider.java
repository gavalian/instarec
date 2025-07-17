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

    public double getM(int pid) {
      switch (pid) {
        case 22:
          return 0;
        case 11:
          return 0.000511;
        case -11:
          return 0.000511;
        case 211:
          return 0.13957;
        case -211:
          return 0.13957;
        case 13:
          return 0.10566;
        case -13:
          return 0.10566;
        case 321:
          return 0.49368;
        case -321:
          return 0.49368;
        case 2212:
          return 0.938272;
        case 2112:
          return 0.939565;
        case 45:
          return 1.875;
        default:
          return -1;
      }
    }

    public double square(double a){
      return a*a;
    }

    public void calcExc(double[] pip, double[] pim, double[] el, double[] exc,double beamE){
      
      double elE=Math.sqrt(square(el[1])+square(getM(11)));
      double pipE=Math.sqrt(square(pip[1])+square(getM(211)));
      double pimE=Math.sqrt(square(pim[1])+square(getM(-211)));
      double pM=getM(2212);

      double IM = Math.sqrt(square(pipE+pimE)- ( square(pip[6]+pim[6]) + square(pip[7]+pim[7]) + square(pip[8]+pim[8]) ));
      double px_m = -1.0*(el[6]+pip[6]+pim[6]);
      double py_m = -1.0*(el[7]+pip[7]+pim[7]);
      double pz_m = beamE - (el[8]+pip[8]+pim[8]);
      double p_m=Math.sqrt(square(px_m) + square(py_m)+square(pz_m));
      double pxp_m = px_m/p_m;
      double pyp_m = py_m/p_m;
      double E_m = beamE + getM(2212) - (elE+pipE+pimE);
      double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));

      exc[0]=IM;
      exc[1]=Math.sqrt(MM2);
      exc[2]=Math.sqrt(pxp_m*pxp_m + pyp_m*pyp_m);
      exc[3]=(MM2-pM*pM)/(2*pM);
    }

    public void getPart(Bank RECPart, int pindex, double[] part){
      if(pindex!=999){
        double pz = RECPart.getFloat("pz", pindex);
        double px = RECPart.getFloat("px", pindex);
        double py = RECPart.getFloat("py", pindex);
        double p=Math.sqrt(px*px+py*py+pz*pz);
        double Theta = Math.acos(pz / p)*(180/Math.PI);// Math.atan2(Math.sqrt(px*px+py*py),pz);
        double Phi = Math.atan2(py, px)*(180/Math.PI);
        int pid=RECPart.getInt("pid", pindex);
        int status=RECPart.getInt("status", pindex);
        int charge=RECPart.getInt("charge", pindex);
        double vz = RECPart.getFloat("vz", pindex);
        part[0]=pid;
        part[1]=p;
        part[2]=Theta;
        part[3]=Phi;
        part[4]=status;
        part[5]=charge;
        part[6]=px;
        part[7]=py;
        part[8]=pz;
        part[9]=vz;
      } else{
        part[0]=999;
        part[1]=999;
        part[2]=999;
        part[3]=999;
        part[4]=999;
        part[5]=999;
        part[6]=999;
        part[7]=999;
        part[8]=999;
        part[9]=999;
      }
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

    public Boolean passFid(double[] Ls, double cut){
      if(Ls[0]>cut && Ls[1]>cut && Ls[2]>cut){
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
    
    public int process(String file, int limTotNegEvs,double threshold, String endName,String output,Boolean clearOutFile){

      //remove existing output files so we can overwrite output
      if(clearOutFile){
        delete_existing_output(output+".csv");
        //if did this per sector
        //for(int i=1;i<14;i++){
        //  if(i<7){delete_existing_output(output+"_sector"+String.valueOf(i)+".csv");}
        //}
      }
      

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


      double[] electron = new double[10];
      double[] elEs=new double[7];
      double[] elLs=new double[9];
      double[] photon = new double[10];
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
        
        StringBuilder csvLineBuilderOut = new StringBuilder();;
        Boolean hasRad=false;

        //same rows in part bank and pred_part bank
        for(int row=0;row<pred_part.getRows();row++){
          short charge = pred_part.getShort(4,row);
          short sector = pred_part.getShort(3,row);
          int pindex=part.getShort(0,row);
          float resp=(float)pred_part.getDouble(2,row);

          double totHTCC=pred_part.getDouble(30,row)+pred_part.getDouble(31,row)+pred_part.getDouble(32,row);

          StringBuilder csvLineBuilder= makeOutString(pred_part, pred_HTCC, pred_ECAL, row, sector);
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
              //remove el rec pid with && recpid!=11
              if(rowphoton!=pindex && recpidph==22 && electron[1]>0.01){

                
                getCalInfo(reccal, pindex, electron[1], elEs, elLs);
                getCalInfo(reccal, rowphoton, photon[1], phEs, phLs);
        
                //make sure part is good electron candidate
                //fids && passFid(elLs, 9)
                //take those not in fid region?
                if(elEs[3]>0.01  && totHTCC!=0 && Math.abs(dTheta)<2.5 && !passFid(elLs, 9)){// elEs[0]>0.06 elEs[5] > (0.15 - elEs[4]) && elEs[3]>0.15

                  
                  double difLU=elLs[0]-phLs[0];
                  double difLV=elLs[1]-phLs[1];
                  double difLW=elLs[2]-phLs[2];
                  double dist=Math.sqrt(difLU*difLU+difLV*difLV+difLW*difLW);

                
                  //region with less bg from non rad photons
                  if( dPhi>-30 && dPhi<30){ //dist>15 &&
                    
                    if(Math.abs(dTheta)<0.5 ){
                      csvLineBuilderOut=csvLineBuilder;
                      hasRad=true;
                      if(resp>threshold){
                        countGoodInsta++;
                      } else {
                        countBadInsta++;
                      }
                      countRad++;
                    } //good dTheta 
                  }//good dPhi
                }//has good electron

              }// has good photon & electron
            }//loop over rec particles to find photon
          }//find matched negative in FD
        }//loop over predicted particles

        if(hasRad==true){
          try{
            FileWriter writer = new FileWriter(output+".csv",true);
            PrintWriter pw = new PrintWriter(writer);
            pw.println(csvLineBuilderOut.toString());
            pw.close();

            //don't count if set limTotNegEvs=-1
            if(limTotNegEvs!=-1){
              count++;
            }
          } catch (IOException e){
            e.printStackTrace();
          }
        }// if has rad, write it out
        
      }//while read

      float eff=(float)countGoodInsta/countRad;

      System.out.printf("\n\nNumber of missIDed electrons with radiated photon %d , recovered %d, not recovered %d\n",countRad,countGoodInsta,countBadInsta);
      System.out.printf("Recovery Efficiency %f\n",eff);

      return countRad;

    }

    public StringBuilder makeOutString(Leaf part,Leaf pred_HTCC,Leaf pred_ECAL,int row, int sector) {
      StringBuilder csvLineBuilder = new StringBuilder();

      float[] pred_strips = new float[9];
      float[] energy = new float[3];
      float sumE = 0;
      for (int j = 0; j < 3; j++) {
        energy[j] = 0;
      }
      float[] DUs = new float[9];
      for (int j = 0; j < pred_ECAL.getRows(); j++) {
        // match cluster to track
        if (pred_ECAL.getShort(0, j) == row) {
          if (pred_ECAL.getInt(2, j) != 0) {
            pred_strips[pred_ECAL.getInt(2, j) - 1] = (float) pred_ECAL.getDouble(3, j);
            DUs[pred_ECAL.getInt(2, j) - 1] = pred_ECAL.getInt(5, j);
            sumE += pred_ECAL.getDouble(4, j);
            // sum energy in PCAL, ECIN, ECOUT
            if (pred_ECAL.getInt(2, j) < 4) {
              energy[0] += pred_ECAL.getDouble(4, j);
            } else if (pred_ECAL.getInt(2, j) >= 4 && pred_ECAL.getInt(2, j) < 7) {
              energy[1] += pred_ECAL.getDouble(4, j);
            } else if (pred_ECAL.getInt(2, j) >= 7) {
              energy[2] += pred_ECAL.getDouble(4, j);
            }
          }
        }
      }

      csvLineBuilder.append(String.format("%.6f,%.6f,%.6f,", energy[0], energy[1], energy[2]));

      for (int j = 0; j < 9; j++) {
        csvLineBuilder.append(String.format("%.6f,", pred_strips[j]));
      }
      for (int j = 0; j < 9; j++) {
        csvLineBuilder.append(String.format("%.6f,", DUs[j]));
      }

      // add wires
      for (int j = 12; j < 18; j++) {
        csvLineBuilder.append(String.format("%.6f,", part.getDouble(j, row)));
      }
      // add HTCC, has one row per sector
      float sumHTCC = 0;
      for (int j = 0; j < 8; j++) {
        csvLineBuilder.append(String.format("%.6f,", pred_HTCC.getDouble(j + 1, sector - 1)));
        sumHTCC += pred_HTCC.getDouble(j + 1, sector - 1);
      }

      int sectorBefore = sector - 1;
      int sectorAfter = sector + 1;
      if (sector == 1) {
        sectorBefore = 6;
      } else if (sector == 6) {
        sectorAfter = 1;
      }

      for (int j = 0; j < 8; j++) {
        csvLineBuilder.append(String.format("%.6f,", pred_HTCC.getDouble(j + 1, sectorBefore - 1)));
        sumHTCC += pred_HTCC.getDouble(j + 1, sectorBefore - 1);
      }

      for (int j = 0; j < 8; j++) {
        csvLineBuilder.append(String.format("%.6f,", pred_HTCC.getDouble(j + 1, sectorAfter - 1)));
        sumHTCC += pred_HTCC.getDouble(j + 1, sectorAfter - 1);
      }

      return csvLineBuilder;
    }

    public int negativesFromTwoPion(String file, int limTotNegEvs,String endName,String output,double beamE, Boolean clearOutFile){

      //remove existing output files so we can overwrite output
      if(clearOutFile){
        delete_existing_output(output+".csv");
        //if did this per sector
        //for(int i=1;i<14;i++){
        //  if(i<7){delete_existing_output(output+"_sector"+String.valueOf(i)+".csv");}
        //}
      }

      int count = 0, countRad=0, countGoodInsta=0, countBadInsta=0;

      HipoReader r = new HipoReader(file);
      Event ev = new Event();
  
      Leaf part = new Leaf(32, 99, "i", 1200);
      Leaf pred_part = new Leaf(32, 3, "i", 1200);
      Leaf pred_ECAL = new Leaf(32, 1, "i", 1200);
      Leaf pred_HTCC = new Leaf(32, 98, "i", 1200);

      Bank recpart = r.getBank("REC::Particle");
      Bank reccal = r.getBank("REC::Calorimeter");

      double[] exc = new double[4];
      double[] el = new double[10];
      double[] elEs=new double[7];
      double[] elLs=new double[9];
      double[] pim = new double[10];
      double[] pimEs=new double[7];
      double[] pimLs=new double[9];
      double[] pip = new double[10];
      double[] pipEs=new double[7];
      double[] pipLs=new double[9];

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

        Boolean writeEvent=false;
        int rowPim=0, SectorPIm=0;
        
        //same rows in part bank and pred_part bank
        for(int row=0;row<pred_part.getRows();row++){
          for(int row2=0;row2<pred_part.getRows();row2++){
            for(int row3=0;row3<pred_part.getRows();row3++){

              int pindex=part.getShort(0,row);
              int pindexpim=part.getShort(0,row2);
              int pindexpip=part.getShort(0,row3);

              if(pindex!=999 && pindexpip!=999 && pindexpim!=999){
                getPart(recpart, pindex, el);
                getCalInfo(reccal, pindex, el[1], elEs, elLs);
                
                getPart(recpart, pindexpim, pim);
                getCalInfo(reccal, pindexpim, pim[1], pimEs, pimLs);

                
                getPart(recpart, pindexpip, pip);
                getCalInfo(reccal, pindexpip, pip[1], pipEs, pipLs);

                int hasPip=0,hasPim=0,hasElCandi=0;

                if(pip[5]==1 && pip[0]==211){
                  hasPip=1;
                }
                if(pim[5]==-1  ){ //&& pim[0]==-211
                  hasPim=1;
                }

                //require good electron so low probability of pim- being an electron
                if(el[5]==-1 && passFid(elLs, 14) && el[0]==11){ //&& matchel!=-1 el[11] is track chi^2 && el[11]<350 && Math.abs(el[13])<20 && el[12]==6
                  hasElCandi=1;
                }

                if(hasElCandi==1 && hasPim==1 && hasPip==1 && row!=row2){

                  calcExc(pip,pim, el, exc,beamE);
                  //System.out.printf("MM %f\n",exc[1]);

                  if(exc[1]>0.8 && exc[1]<1.1 && el[1]>2){ //exc[0]>0.6 && exc[0]<0.9
                    writeEvent=true;
                    rowPim=row2;
                  }//cuts on exclusivity
                }//required candidates

              }//good pindices
    
            }//loop over part 3
          }//loop over part 2
        }//loop over part1
        if(writeEvent){
          try{
            StringBuilder csvLineBuilder= makeOutString(pred_part, pred_HTCC, pred_ECAL, rowPim, pred_part.getShort(3,rowPim));
            csvLineBuilder.append("1,0");
            FileWriter writer = new FileWriter(output+".csv",true);
            PrintWriter pw = new PrintWriter(writer);
            pw.println(csvLineBuilder.toString());
            pw.close();
            //don't count if set limTotNegEvs=-1
            if(limTotNegEvs!=-1){
              count++;
            }
          } catch (IOException e){
            e.printStackTrace();
          }
        }
      }//loop over events
      System.out.printf("\n\nWrote %d Events\n",count);
      return count;
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

    public Integer copyOverNegatives(int numLines, String inputPath, String outputPath) {
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
            return count;
        } catch (IOException e) {
            System.err.println("Error processing files: " + e.getMessage());
            return -1;
        }
    }

    public Integer copyOverPositives(int numLines, String inputPath, String outputPath) {
      try (
          BufferedReader reader = new BufferedReader(new FileReader(inputPath+".csv"));
          BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath+".csv", true)) // append = true
      ) {
          String line;
          int count = 0;

          while ((line = reader.readLine()) != null && count < numLines) {
              if (line.trim().endsWith("0,1")) {
                  writer.write(line);
                  writer.newLine();
                  count++;
              }
          }
          System.out.println("Appended " + count + " positives to: " + outputPath);
          return count;
      } catch (IOException e) {
          System.err.println("Error processing files: " + e.getMessage());
          return -1;
      }
  }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar -in wvalid.h5 -inneg wtrainNeg.h5 -o training_data/ElPIDTrainRadPhotons -nr training_data/ElPIDTrain
    //read plots in j4shell with eg TwigStudio.browser("plots/RadPhoton0.twig");
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting radiated photon data provider ");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.addRequired("-inneg", "input name");
      p.addRequired("-o", "output file name");
      p.addRequired("-nr", "usual (not rad) training file name");
      p.parse(args);

      String endName="_noFid"; // used to change output path of plots (eg adding _NoFiducialCuts)
        
      RadPhotonDataProvider dp = new RadPhotonDataProvider();
      int nNegatives=dp.negativesFromTwoPion(p.getOption("-inneg").stringValue(),-1,endName,p.getOption("-o").stringValue(),10.6,true);
      int nPositives=dp.process(p.getOption("-in").stringValue(),nNegatives,0.05,endName,p.getOption("-o").stringValue(),false);
      //int nNegatives=dp.copyOverNegatives(nPositives,p.getOption("-nr").stringValue(),p.getOption("-o").stringValue());
        
    }
}
