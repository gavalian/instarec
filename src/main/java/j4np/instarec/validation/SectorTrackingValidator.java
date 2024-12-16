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
import twig.data.StatNumber;
import twig.data.TDirectory;
import twig.graphics.TGCanvas;

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
public class SectorTrackingValidator {
    
    public SectorTrackingValidator(){
    }

    public int findRECTrackSector(Bank trackbank,int pindex){

      for(int row=0;row<trackbank.getRows();row++){
        int sector = trackbank.getInt("sector",row);
        int pind = trackbank.getInt("pindex",row);
        if(pind==pindex){
          return sector;
        }
      }

      return -1;

    }

    public Boolean has6SL(Bank trackbank, int pindex){
      if( (trackbank.getInt("status", pindex) & 0b101010101010)!=0){
        return false;
      } else {
        return true;
      }
    }
    
    public void process(String file, int limEvs, String endName, double threshold, Boolean sixSLOnly){

      double[] InstaSecDif= new double[4];
      double[] RecSecDif= new double[4];
      double[] TP= new double[4];
      double[] FN= new double[4];
      for (int l=0; l<4;l++){
        InstaSecDif[l]=0;
        RecSecDif[l]=0;
        TP[l]=0;
        FN[l]=0;
      }

      //output training sample size per sector
      int count = 0;

      HipoReader r = new HipoReader(file);
      Event ev = new Event();

      Leaf part = new Leaf(32200, 99, "i", 1200);
      Leaf pred_part = new Leaf(32200, 1, "i", 1200);
      Bank recpart = r.getBank("REC::Particle");
      Bank rectrack = r.getBank("REC::Track");

      //don't count if set limEvs=-1
      if(limEvs==-1){count=-2;}
      
      while(r.hasNext() && count<limEvs){
        //don't count if set limEvs=-1
        if(limEvs!=-1){
          count++;
        }

        double[] InstaSecDif_pEv= new double[4];
        double[] RecSecDif_pEv= new double[4];
        for (int l=0; l<4;l++){
          InstaSecDif_pEv[l]=0;
          RecSecDif_pEv[l]=0;
        }

        r.nextEvent(ev);
        ev.read(part,32200,99);
        ev.read(pred_part,32200,1);
        ev.read(recpart);
        ev.read(rectrack);
        for(int row=0;row<pred_part.getRows();row++){
          short charge = pred_part.getShort(4,row);
          short sector = pred_part.getShort(3,row);
          double resp=pred_part.getDouble(2,row);
          //ElPIDWorker set resp to 0.5 for pos tracks
          if(charge==1){resp=0;}
          //check that we have non electron FD track
          if(sector>0 && resp<threshold){
            for (int row2 = row + 1; row2 < pred_part.getRows(); row2++) {
              short charge2 = pred_part.getShort(4,row2);
              short sector2 = pred_part.getShort(3,row2);
              double resp2=pred_part.getDouble(2,row2);
              //ElPIDWorker set resp to 0.5 for pos tracks
              if(charge2==1){resp2=0;}
              //check that we have non electron FD track
              if(sector2>0 && resp2<threshold){
                int sDif=(Math.abs(sector-sector2));
                if(sDif==4){sDif=2;};
                if(sDif==5){sDif=1;};
                InstaSecDif_pEv[sDif]++;
              } //check that second track is not an electron
            }//find other predicted track in event
          }//check that it's not an electron in FD
        }//loop over predicted tracks

        for(int row=0;row<recpart.getRows();row++){
          int charge = recpart.getInt("charge",row);
          int sector=findRECTrackSector(rectrack,row);
          short status=recpart.getShort("status",row);
          int pid=recpart.getInt("pid",row);
          float vz = recpart.getFloat("vz", row);
          Boolean rightNSL=true;
          if(sixSLOnly){
            rightNSL=has6SL(rectrack,row);
          }
          //check that we have good non electron FD track
          if(sector>0 && pid!=11 && charge!=0 && Math.abs(vz)<20 && Math.abs(status)>=2000 && Math.abs(status)<4000 && rightNSL && rectrack.getFloat("chi2", row)<350){
            for (int row2 = row + 1; row2 < recpart.getRows(); row2++) {
              int charge2 = recpart.getInt("charge",row2);
              int sector2 = findRECTrackSector(rectrack,row2);
              short status2=recpart.getShort("status",row2);
              int pid2=recpart.getInt("pid",row2);
              float vz2 = recpart.getFloat("vz", row2);
              Boolean rightNSL2=true;
              if(sixSLOnly){
                rightNSL2=has6SL(rectrack,row2);
              }
              //check that we have good non electron FD track
              if(sector2>0 && pid2!=11 && charge2!=0 && Math.abs(vz2)<20 && Math.abs(status2)>=2000 && Math.abs(status2)<4000 && rightNSL2 && rectrack.getFloat("chi2", row2)<350){
                int sDif=(Math.abs(sector-sector2));
                if(sDif==4){sDif=2;};
                if(sDif==5){sDif=1;};
                RecSecDif_pEv[sDif]++;
              } //check that second track is not an electron
            }//find other REC track in event
          }//check that it's not an electron in FD
        }//loop over REC tracks

        for (int l=0; l<4;l++){
          if(RecSecDif_pEv[l]>0){
            RecSecDif[l]+=RecSecDif_pEv[l];

            if(InstaSecDif_pEv[l]>=0){
              InstaSecDif[l]+=InstaSecDif_pEv[l];
            }
            
            if(InstaSecDif_pEv[l]>=RecSecDif_pEv[l]){
              TP[l]+=InstaSecDif_pEv[l];
            } else{
              FN[l]+=InstaSecDif_pEv[l];
            }
          }
        }

      }//loop over events

      GraphErrors gEff = new GraphErrors();
      gEff.attr().setMarkerColor(4);
      gEff.attr().setMarkerSize(15);
      gEff.attr().setLineWidth(3);
      gEff.attr().setTitle("Efficiency vs Sector Difference");
      gEff.attr().setTitleX("Sector Difference");
      gEff.attr().setTitleY("Efficiency");

      for (int l=0; l<4;l++){

        /*StatNumber Eff = new StatNumber(TP[l], Math.sqrt(TP[l]));
        StatNumber FNs = new StatNumber(FN[l], Math.sqrt(FN[l]));
        StatNumber denomEff = new StatNumber(TP[l], Math.sqrt(TP[l]));
        denomEff.add(FNs);
        Eff.divide(denomEff);

        System.out.printf("Sector Difference %d: Efficiency %f +/- %f\n",l,Eff.number(),Eff.error());
        
        gEff.addPoint(l, Eff.number(), 0,Eff.error());*/

        double Eff=TP[l]/(TP[l]+FN[l]);
        double EffErr=Math.sqrt(((TP[l]+1)/(TP[l]+FN[l]+2))*((TP[l]+2)/(TP[l]+FN[l]+3))-((TP[l]+1)/(TP[l]+FN[l]+2))*((TP[l]+1)/(TP[l]+FN[l]+2)));
        System.out.printf("Sector Difference %d: Efficiency %f +/- %f\n",l,Eff,EffErr);
        
        gEff.addPoint(l, Eff, 0,EffErr);
      }

      TDirectory.export("plots/Tracking" + endName + ".twig", "/ai/validation/SectorDifference", gEff);

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
    //read plots in j4shell with eg TwigStudio.browser("plots/Tracking.twig");
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting sector difference tracking validator ");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.parse(args);

      String endName=""; // used to change output path of plots (eg adding _NoFiducialCuts)

        
      SectorTrackingValidator dp = new SectorTrackingValidator();
      dp.process(p.getOption("-in").stringValue(),1500000,endName,0.075,true);
  
      
      
        
    }
}
