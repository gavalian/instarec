/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.validation;

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
import java.io.IOException;

//for priting
import j4np.hipo5.data.Bank;

/**
 *
 * @author tyson
 */
public class RadPhotonValidator {
    
    public RadPhotonValidator(){
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
    
    public void process(String file, int limTotNegEvs,double threshold, String endName){

      //output training sample size per sector
      int count = 0, countRad=0, countGoodInsta=0, countBadInsta=0;

      HipoReader r = new HipoReader(file);
      Event ev = new Event();
  
      Leaf part = new Leaf(32200, 99, "i", 1200);
      Leaf pred_part = new Leaf(32200, 1, "i", 1200);

      Bank recpart = r.getBank("REC::Particle");
      Bank reccal = r.getBank("REC::Calorimeter");

      H2F hThetaPhi = new H2F("dPhi vs dTheta",50,-1.0,1.0,50,-15,5);
      hThetaPhi.attr().setTitleX("dTheta [Degrees]");
      hThetaPhi.attr().setTitleY("dPhi [Degrees]");

      H2F hSfP = new H2F("Sampling Fraction vs P",70,1,8,50,0,0.5);
      hSfP.attr().setTitleX("P [GeV]");
      hSfP.attr().setTitleY("Sampling Fraction");

      H2F hEP = new H2F("PCAL Energy vs P",70,1,8,50,0,0.15);
      hEP.attr().setTitleX("P [GeV]");
      hEP.attr().setTitleY("PCAL Energy Deposition [GeV]");

      H2F hDistdTheta = new H2F("Distance between Photon and Electron vs dTheta",50,-1.0,1.0,50,0,100);
      hDistdTheta.attr().setTitleX("dTheta [Degrees]");
      hDistdTheta.attr().setTitleY("Distance [cm]");

      H1F hDist = new H1F("Distance between Photon and Electron",50,0,100);
      hDist.attr().setTitleX("Distance [cm]");
      hDist.attr().setLineColor(2);
      hDist.attr().setFillColor(2);
      hDist.attr().setLineWidth(3);

      H1F hTheta = new H1F("All",50,-5,5);
      hTheta.attr().setTitleX("dTheta [Degrees]");
      hTheta.attr().setLineColor(2);
      hTheta.attr().setFillColor(2);
      hTheta.attr().setLineWidth(3);
      
      H1F hTheta2 = new H1F("With Online e-",50,-5,5);
      hTheta2.attr().setTitleX("dTheta [Degrees]");
      hTheta2.attr().setLineColor(5);
      hTheta2.attr().setLineWidth(3);

      H1F hP_Insta = new H1F("Offline ! e- PID & Online e- PID",50,0.,10);
      hP_Insta.attr().setLineColor(2);
      hP_Insta.attr().setFillColor(2);
      hP_Insta.attr().setLineWidth(3);
      hP_Insta.attr().setTitleX("Momentum [GeV]");
      H1F hP_InstaBad = new H1F("Offline ! e- PID & Online ! e- PID",50,0.,10);
      hP_InstaBad.attr().setLineColor(3);
      hP_InstaBad.attr().setLineWidth(3);
      hP_InstaBad.attr().setTitleX("Momentum [GeV]");
      H1F hP_REC = new H1F("Offline ! e- PID",50,0.,10);
      hP_REC.attr().setLineColor(5);
      hP_REC.attr().setLineWidth(3);
      hP_REC.attr().setTitleX("Momentum [GeV]");



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
        ev.read(part,32200,99);
        ev.read(pred_part,32200,1);
        ev.read(recpart);
        ev.read(reccal);
        

        //same rows in part bank and pred_part bank
        for(int row=0;row<pred_part.getRows();row++){
          short charge = pred_part.getShort(4,row);
          short sector = pred_part.getShort(3,row);
          int pindex=part.getShort(0,row);
          float resp=(float)pred_part.getDouble(2,row);

          double totHTCC=pred_part.getDouble(30,row)+pred_part.getDouble(31,row)+pred_part.getDouble(32,row);

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
              if(rowphoton!=pindex && recpid!=11 && recpidph==22 ){

                
                getCalInfo(reccal, pindex, electron[1], elEs, elLs);
                getCalInfo(reccal, rowphoton, photon[1], phEs, phLs);
  
                hSfP.fill(electron[1],elEs[3]);
                hEP.fill(electron[1],elEs[0]);
        
                //make sure part is good electron candidate
                if(elEs[0]>0.06 && elEs[3]>0.15 && totHTCC!=0 && !passFid(elLs)){// elEs[5] > (0.15 - elEs[4]) 

                  hThetaPhi.fill(dTheta, dPhi);
                  
                  double difLU=elLs[0]-phLs[0];
                  double difLV=elLs[1]-phLs[1];
                  double difLW=elLs[2]-phLs[2];
                  double dist=Math.sqrt(difLU*difLU+difLV*difLV+difLW*difLW);

                  hDist.fill(dist);
                  hDistdTheta.fill(dTheta,dist);
                
                  //region with less bg from non rad photons
                  if(dist>30 && dPhi>-30 && dPhi<30){
                    //don't count if set limTotNegEvs=-1
                    if(limTotNegEvs!=-1){
                      count++;
                    }
                    hTheta.fill(dTheta);

                    if(resp>threshold){
                      hTheta2.fill(dTheta);
                      
                    }
                    if(Math.abs(dTheta)<0.5 ){
                      countRad++;
                      hP_REC.fill(electron[1]);
                      if(resp>threshold){
                        countGoodInsta++;
                        hP_Insta.fill(electron[1]);
                      } else {
                        countBadInsta++;
                        hP_InstaBad.fill(electron[1]);
                      }

                    }  
                  }
                }

              }
            }
          }
        }
      }

      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/ThetaDif", hThetaPhi);
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/ThetaDif", hTheta);
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/ThetaDif", hTheta2);
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/ValidVars", hDist);
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/ValidVars", hDistdTheta); 
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/ValidVars", hSfP);
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/ValidVars", hEP);
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/elP", hP_REC);
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/elP", hP_Insta);
      TDirectory.export("plots/RadPhoton0" + endName + ".twig", "/ai/validation/elP", hP_InstaBad);

      float eff=(float)countGoodInsta/countRad;

      System.out.printf("\n\nNumber of missIDed electrons with radiated photon %d , recovered %d, not recovered %d\n",countRad,countGoodInsta,countBadInsta);
      System.out.printf("Recovery Efficiency %f\n",eff);

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
    //read plots in j4shell with eg TwigStudio.browser("plots/RadPhoton0.twig");
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting radiated photon validator ");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.parse(args);

      String endName="_antiFid"; // used to change output path of plots (eg adding _NoFiducialCuts)
        
      RadPhotonValidator dp = new RadPhotonValidator();
      dp.process(p.getOption("-in").stringValue(),-1,0.075,endName);
  
      
      
        
    }
}
