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

    public static int[] convertL1Trigger(long bits) {
      int[] trigger = new int[32];
  
      // System.out.printf("%X - %X\n", bits,bits&0xF);
      for (int i = 0; i < trigger.length; i++) {
        trigger[i] = 0;
        if (((bits >> i) & (1L)) != 0L)
          trigger[i] = 1;
        // System.out.println(Arrays.toString(trigger));
      }
      return trigger;
    }

    public int hasTriggerMesonEx(Bank triggerBank){
      long bits = triggerBank.getLong("trigger", 0);
      int[] L1trigger = convertL1Trigger(bits);
      return L1trigger[25];
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

      //Leaf pred_part = new Leaf(32200, 1, "i", 1200);
      Bank recpart = r.getBank("REC::Particle");
      Bank rectrack = r.getBank("REC::Track");
      Bank triggerbank = r.getBank("RUN::config");

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
        Leaf pred_part = ev.readLeaf(1,12,32200,1);
        Leaf pred_track = ev.readLeaf(1,12,32000,1);
        //ev.read(pred_part,32200,1);
        ev.read(recpart);
        ev.read(rectrack);
        ev.read(triggerbank);

        int hasTrig=hasTriggerMesonEx(triggerbank);

        //System.out.println("\n\nnew event: ");
        //pred_track.print();
        //System.out.println("\nparticles: \n");
        //pred_part.print();
        //triggerbank.show();


        if(hasTrig==1 && recpart.getRows()>0){
          for(int row=0;row<pred_part.getRows();row++){
            short charge = pred_part.getShort(4,row);
            short sector = pred_part.getShort(3,row);
            double resp=pred_part.getDouble(2,row);
            short pindex = pred_part.getShort(0,row);
            short goodtrack = pred_track.getShort(0,pindex);

            
            //ElPIDWorker set resp to 0.5 for pos tracks
            if(charge==1){resp=0;}
            //check that we have non electron FD track
            if(sector>0 && goodtrack>0){ // && resp<threshold){
              for (int row2 = row + 1; row2 < pred_part.getRows(); row2++) {
                short charge2 = pred_part.getShort(4,row2);
                short sector2 = pred_part.getShort(3,row2);
                double resp2=pred_part.getDouble(2,row2);
                short pindex2 = pred_part.getShort(0,row2);
                short goodtrack2 = pred_track.getShort(0,pindex2);
                //System.out.printf("pindex %d pindex2 %d gc %d gc2 %d \n",pindex,pindex2,goodtrack,goodtrack2);
                //ElPIDWorker set resp to 0.5 for pos tracks
                if(charge2==1){resp2=0;}
                //check that we have non electron FD track
                if(sector2>0 && goodtrack2>0 ){ // && resp2<threshold){
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
            if(sector>0  && charge!=0 && Math.abs(vz)<20 && Math.abs(status)>=2000 && Math.abs(status)<4000 && rightNSL && rectrack.getFloat("chi2", row)<350){ //&& pid!=11
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
                if(sector2>0 && charge2!=0 && Math.abs(vz2)<20 && Math.abs(status2)>=2000 && Math.abs(status2)<4000 && rightNSL2 && rectrack.getFloat("chi2", row2)<350){ // && pid2!=11
                  int sDif=(Math.abs(sector-sector2));
                  if(sDif==4){sDif=2;};
                  if(sDif==5){sDif=1;};
                  RecSecDif_pEv[sDif]++;
                } //check that second track is not an electron
              }//find other REC track in event
            }//check that it's not an electron in FD
          }//loop over REC tracks
        } //check has trig

        for (int l=0; l<4;l++){
          RecSecDif[l]+=RecSecDif_pEv[l];
          InstaSecDif[l]+=InstaSecDif_pEv[l];
        }

      }//loop over events

      GraphErrors gEff = new GraphErrors();
      gEff.attr().setMarkerColor(4);
      gEff.attr().setMarkerSize(15);
      gEff.attr().setLineWidth(3);
      gEff.attr().setTitle("Efficiency vs Sector Difference");
      gEff.attr().setTitleX("Sector Difference");
      gEff.attr().setTitleY("Efficiency");

      GraphErrors gNTracks = new GraphErrors();
      gNTracks.attr().setMarkerColor(2);
      gNTracks.attr().setMarkerSize(15);
      gNTracks.attr().setLineWidth(3);
      gNTracks.attr().setTitle("Online Tracking");
      gNTracks.attr().setTitleX("Sector Difference");
      gNTracks.attr().setTitleY("Number of Tracks");

      GraphErrors gNTracks_rec = new GraphErrors();
      gNTracks_rec.attr().setMarkerColor(5);
      gNTracks_rec.attr().setMarkerSize(15);
      gNTracks_rec.attr().setLineWidth(3);
      gNTracks_rec.attr().setTitle("Offline Tracking");
      gNTracks_rec.attr().setTitleX("Sector Difference");
      gNTracks_rec.attr().setTitleY("Number of Tracks");

      for (int l=0; l<4;l++){

        double Eff=InstaSecDif[l]/(RecSecDif[l]);
        double EffErr=Math.sqrt(((RecSecDif[l]+1)/(RecSecDif[l]+2))*((RecSecDif[l]+2)/(RecSecDif[l]+3))-((RecSecDif[l]+1)/(RecSecDif[l]+2))*((RecSecDif[l]+1)/(RecSecDif[l]+2)));
        System.out.printf("Sector Difference %d: Efficiency %f +/- %f\n",l,Eff,EffErr);
        
        gEff.addPoint(l, Eff, 0,EffErr);
        gNTracks.addPoint(l,InstaSecDif[l],0,Math.sqrt(InstaSecDif[l]));
        gNTracks_rec.addPoint(l,RecSecDif[l],0,Math.sqrt(RecSecDif[l]));
      }

      TDirectory.export("plots/Tracking" + endName + ".twig", "/ai/validation/SectorDifference_efficiency", gEff);
      TDirectory.export("plots/Tracking" + endName + ".twig", "/ai/validation/SectorDifference_NTracks", gNTracks);
      TDirectory.export("plots/Tracking" + endName + ".twig", "/ai/validation/SectorDifference_NTracksRec", gNTracks_rec);

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
      dp.process(p.getOption("-in").stringValue(),100000,endName,0.075,false);
  
      
      
        
    }
}
