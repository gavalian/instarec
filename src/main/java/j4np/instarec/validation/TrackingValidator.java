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
public class TrackingValidator {
    
    public TrackingValidator(){
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

    public int fillMatchPredTrack(Leaf track, double[] part, double[] pred_part, int desired_charge, double lim_p) {

      pred_part[0]=999;
      short best_pindex=999;
      double best_resp=999;
      for (short i = 0; i < track.getRows(); i++) {
        short charge = track.getShort(3,i);
        double px = track.getDouble(5,i);
        double py = track.getDouble(6,i);
        double pz = track.getDouble(7,i);
        double p=Math.sqrt(px*px+py*py+pz*pz);
        double theta = (Math.acos(pz / p)*(180/Math.PI));// Math.atan2(Math.sqrt(px*px+py*py),pz);
        double phi = (Math.atan2(py, px)*(180/Math.PI));
        double resp=Math.abs(part[1]-p);
        if(resp<best_resp && resp<(lim_p*part[1]) && desired_charge==charge){ 
          best_pindex=i;
          best_resp=resp;
          pred_part[0]=resp;
          pred_part[1]=p;
          pred_part[2]=theta;
          pred_part[3]=phi;
        }
      }

      return best_pindex;
    }

    public StatNumber[] calcMets(float TP, float FP, float FN){
      //StatNumber Pur = new StatNumber(TP, Math.sqrt(TP));
      StatNumber Eff = new StatNumber(TP, Math.sqrt(TP));
      //StatNumber FPs = new StatNumber(FP, Math.sqrt(FP));
      StatNumber FNs = new StatNumber(FN, Math.sqrt(FN));
      //StatNumber denomPur = new StatNumber(TP, Math.sqrt(TP));
      StatNumber denomEff = new StatNumber(TP, Math.sqrt(TP));
      //denomPur.add(FPs);
      denomEff.add(FNs);
      //Pur.divide(denomPur);
      Eff.divide(denomEff);
      StatNumber[] mets = new StatNumber[2];
      //mets[0] = Pur;
      mets[1] = Eff;
      return mets;
    }

    public Boolean has6SL(Bank trackbank, int pindex){
      if( (trackbank.getInt("status", pindex) & 0b101010101010)!=0){
        return false;
      } else {
        return true;
      }
    }

    public void makeGraphs(int nbins, float[] bincentre, float[] TP, float[] FP, float[] FN,String name, String title, String units, String ChargeSt, String PIDSt, String endName){
      GraphErrors gEff = new GraphErrors();
      gEff.attr().setMarkerColor(2);
      gEff.attr().setMarkerSize(10);
      gEff.attr().setTitle("Efficiency vs "+title);
      gEff.attr().setTitleX(title+" "+units);
      gEff.attr().setTitleY("Metrics");
      GraphErrors gPur = new GraphErrors();
      gPur.attr().setMarkerColor(5);
      gPur.attr().setMarkerSize(10);
      gPur.attr().setTitle("Purity vs "+title);
      gPur.attr().setTitleX(title+" "+units);
      gPur.attr().setTitleY("Metrics");

      for(int i=0;i<nbins;i++){
        StatNumber mets[]=calcMets(TP[i], FP[i], FN[i]);
        //double Pur = mets[0].number();
        double Eff = mets[1].number();
        //double PurErr = mets[0].error();
        double EffErr = mets[1].error();
        //gPur.addPoint(bincentre[i], Pur, 0, PurErr);
        gEff.addPoint(bincentre[i], Eff, 0, EffErr);
      }

      //TDirectory.export("plots/Tracking" + endName + ".twig", "/ai/validation/Purity_"+name, gPur);
      TDirectory.export("plots/Tracking" + endName + ".twig", "/ai/validation/Efficiency_"+name+" ( "+ChargeSt+" Charge, PID "+PIDSt+" )", gEff);

    }
    
    public void process(String file, int limEvs, float[] pbins, float[] thetabins, float[] phibins,String endName, int desired_charge, int desired_pid, double lim_p, Boolean sixSLOnly, Boolean print){

      float TP=0,FP=0,FN=0;

      int npBins=0;
      for(float bin=pbins[1];bin<(pbins[2]-pbins[0]);bin+=pbins[0]){npBins++;}
      float[] bincentre_p=new float[npBins];
      float[] TP_p=new float[npBins];
      float[] FP_p=new float[npBins];
      float[] FN_p=new float[npBins];
      for(int i=0;i<npBins;i++){
        bincentre_p[i]=0;
        TP_p[i]=0;
        FP_p[i]=0;
        FN_p[i]=0;
      }

      int nthetaBins=0;
      for(float bin=thetabins[1];bin<(thetabins[2]-thetabins[0]);bin+=thetabins[0]){nthetaBins++;}
      float[] bincentre_theta=new float[nthetaBins];
      float[] TP_theta=new float[nthetaBins];
      float[] FP_theta=new float[nthetaBins];
      float[] FN_theta=new float[nthetaBins];
      for(int i=0;i<nthetaBins;i++){
        bincentre_theta[i]=0;
        TP_theta[i]=0;
        FP_theta[i]=0;
        FN_theta[i]=0;
      }

      int nphiBins=0;
      for(float bin=phibins[1];bin<(phibins[2]-phibins[0]);bin+=phibins[0]){nphiBins++;}
      float[] bincentre_phi=new float[nphiBins];
      float[] TP_phi=new float[nphiBins];
      float[] FP_phi=new float[nphiBins];
      float[] FN_phi=new float[nphiBins];
      for(int i=0;i<nphiBins;i++){
        bincentre_phi[i]=0;
        TP_phi[i]=0;
        FP_phi[i]=0;
        FN_phi[i]=0;
      }

      //output training sample size per sector
      int count = 0;

      HipoReader r = new HipoReader(file);
      Event ev = new Event();

      /*for priting */
      Leaf trackbank = new Leaf(32000,1,"i",4096);
      Bank recpart = r.getBank("REC::Particle");
      Bank rectrack = r.getBank("REC::Track");

      // pid, P, Theta, Phi
      double[] part=new double[4];
      double[] pred_part=new double[4];

      //don't count if set limEvs=-1
      if(limEvs==-1){count=-2;}
      
      while(r.hasNext() && count<limEvs){

        r.nextEvent(ev);
        ev.read(trackbank,32000,1);
        ev.read(recpart);
        ev.read(rectrack);

        float FN_rec=FN;

        for(int i=0; i<recpart.getRows();i++){
          int charge=recpart.getInt("charge",i);
          short status=recpart.getShort("status",i);
          float vz = recpart.getFloat("vz", i);
          float chi2pid = recpart.getFloat("chi2pid", i);
          getPart(recpart, i, part);

          Boolean rightNSL=true;
          if(sixSLOnly){
            rightNSL=has6SL(rectrack,i);
          }

          //allow to get efficiency for all pid/charge
          if(desired_pid==0){
            //don't want to use 0 pid as bad parts
            if(part[0]!=0){
              part[0]=desired_pid;
              chi2pid=0;
            } else{
              part[0]=999;
            }
          }
          if(desired_charge==0){
            //don't want to use neutral parts as no tracks
            if(charge!=0){
              charge=desired_charge;
            } else{
              charge=999;
            }
          }

          //good track in FD
          if(Math.abs(vz)<20 && rectrack.getFloat("chi2", i)<350 && rightNSL && Math.abs(status)>=2000 && Math.abs(status)<4000){
            //get efficiency for certain pid type or charge
            if(part[0]==desired_pid && Math.abs(chi2pid)<5 && charge==desired_charge){

              //don't count if set limEvs=-1
              if(limEvs!=-1){
                count++;
              }

              int match=fillMatchPredTrack(trackbank, part, pred_part, charge, lim_p);

              if(match!=999){
                TP++;
              } else{
                FN++;

                if(print){
                  System.out.println("\nMissed track in this event:");
                  System.out.printf("pid %d px %f py %f pz %f\n",recpart.getInt("pid",i),recpart.getFloat("px",i),recpart.getFloat("py",i),recpart.getFloat("pz",i));
                  trackbank.print();

                }
              }

              int nbin=0;
              for(float bin=pbins[1];bin<(pbins[2]-pbins[0]);bin+=pbins[0]){
                if(part[1]>=bin && part[1]<(bin+pbins[0])){
                  bincentre_p[nbin]=bin+(pbins[0]/2);
                  if(match!=999){
                    TP_p[nbin]++;
                  } else {
                    FN_p[nbin]++;
                  }
                }
                nbin++;
              }
  
              nbin=0;
              for(float bin=thetabins[1];bin<(thetabins[2]-thetabins[0]);bin+=thetabins[0]){
                if(part[2]>=bin && part[2]<(bin+thetabins[0])){
                  bincentre_theta[nbin]=bin+(thetabins[0]/2);
                  if(match!=999){
                    TP_theta[nbin]++;
                  } else {
                    FN_theta[nbin]++;
                  }
                }
                nbin++;
              }
  
              nbin=0;
              for(float bin=phibins[1];bin<(phibins[2]-phibins[0]);bin+=phibins[0]){
                if(part[3]>=bin && part[3]<(bin+phibins[0])){
                  bincentre_phi[nbin]=bin+(phibins[0]/2);
                  if(match!=999){
                    TP_phi[nbin]++;
                  } else {
                    FN_phi[nbin]++;
                  }
                }
                nbin++;
              }


            }
          }

        }

        /*if(FN>FN_rec && print){
          System.out.println("\nMissed some tracks in this event:");
          recpart.show();
          trackbank.print();
        }*/

      }

      double eff=TP/(TP+FN);
      String ChargeSt="Negative";
      if(desired_charge==1){ChargeSt="Positive";}
      if(desired_charge==0){ChargeSt="All";}

      String PIDSt=String.valueOf(desired_pid);
      if(desired_pid==0){PIDSt="All";}

      System.out.println("\n\nFor "+ChargeSt+" Charge & "+PIDSt+" PID: ");
      System.out.printf("TP %f FN %f \n",TP,FN);
      System.out.printf("eg Efficiency of %f \n",eff);

      makeGraphs(npBins,bincentre_p,TP_p,FP_p,FN_p,"P","Momentum","[GeV]", ChargeSt, PIDSt, endName);
      makeGraphs(nthetaBins,bincentre_theta,TP_theta,FP_theta,FN_theta,"theta","Theta","[Degrees]", ChargeSt, PIDSt, endName);
      makeGraphs(nphiBins,bincentre_phi,TP_phi,FP_phi,FN_phi,"phi","Phi","[Degrees]", ChargeSt, PIDSt, endName);

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
        
      System.out.println("\n\n----- starting tracking validator ");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.parse(args);

      String endName=""; // used to change output path of plots (eg adding _NoFiducialCuts)

      float[] respbins = new float[3];
      respbins[0]=(float)0.01;
      respbins[1]=(float)0.01;
      respbins[2]=(float)0.99;

      float[] pbins = new float[3];
      pbins[0]=(float)1;
      pbins[1]=(float)1;
      pbins[2]=(float)9;

      float[] thetabins = new float[3];
      thetabins[0]=(float)5;
      thetabins[1]=(float)5;
      thetabins[2]=(float)35;

      float[] phibins = new float[3];
      phibins[0]=(float)10;
      phibins[1]=(float)-180;
      phibins[2]=(float)180;

      double lim_p_res=0.2; //percentage
        
      TrackingValidator dp = new TrackingValidator();
      dp.process(p.getOption("-in").stringValue(),150000,pbins,thetabins,phibins,endName,-1,0,lim_p_res,true,false);
      dp.process(p.getOption("-in").stringValue(),150000,pbins,thetabins,phibins,endName,1,0,lim_p_res,true,false);
      dp.process(p.getOption("-in").stringValue(),150000,pbins,thetabins,phibins,endName,-1,11,lim_p_res,true,false);
  
      
      
        
    }
}
