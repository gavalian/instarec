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
import twig.data.GraphErrors;
import twig.data.H1F;
import twig.data.StatNumber;
import twig.data.TDirectory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.BufferedReader;
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
public class ElPIDValidater {
    
    public ElPIDValidater(){
    }

    public int[] checkmets(float resp, double threshold, int pid,double htcc){
      int[] mets=new int[3];
      mets[0]=0;
      mets[1]=0;
      mets[2]=0;
      if(resp>threshold){
        if(pid==11 && htcc!=0){ 
          mets[0]=1;
        } else{
          mets[1]=1;
        }
      } else{
        if(pid==11 && htcc!=0){ 
          mets[2]=1;
        }
      }
      return mets;
    }

    public StatNumber[] calcMets(float TP, float FP, float FN){
      StatNumber Pur = new StatNumber(TP, Math.sqrt(TP));
      StatNumber Eff = new StatNumber(TP, Math.sqrt(TP));
      StatNumber FPs = new StatNumber(FP, Math.sqrt(FP));
      StatNumber FNs = new StatNumber(FN, Math.sqrt(FN));
      StatNumber denomPur = new StatNumber(TP, Math.sqrt(TP));
      StatNumber denomEff = new StatNumber(TP, Math.sqrt(TP));
      denomPur.add(FPs);
      denomEff.add(FNs);
      Pur.divide(denomPur);
      Eff.divide(denomEff);
      StatNumber[] mets = new StatNumber[2];
      mets[0] = Pur;
      mets[1] = Eff;
      return mets;
    }

    public void makeGraphs(int nbins, float effLow,float[] bincentre, float[] TP, float[] FP, float[] FN,String name, String title, String units){
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

      double bestPuratEffLow=0;
      float bestBin=0;


      for(int i=0;i<nbins;i++){
        StatNumber mets[]=calcMets(TP[i], FP[i], FN[i]);
        double Pur = mets[0].number();
        double Eff = mets[1].number();
        double PurErr = mets[0].error();
        double EffErr = mets[1].error();
        gPur.addPoint(bincentre[i], Pur, 0, PurErr);
        gEff.addPoint(bincentre[i], Eff, 0, EffErr);
        if (Eff > effLow) {
          if (Pur > bestPuratEffLow) {
            bestPuratEffLow = Pur;
            bestBin = bincentre[i];
          }
        }
      }

      System.out.format("%n Best Purity at Efficiency above %f: %.3f at %.3f in "+title+" "+units+" \n\n",
        effLow, bestPuratEffLow, bestBin);

      TDirectory.export("plots/ElPID0" + ".twig", "/ai/validation/Purity_"+name, gPur);
      TDirectory.export("plots/ElPID0" + ".twig", "/ai/validation/Efficiency_"+name, gEff);

    }
    
    public void process(String file, int limTotNegEvs,double threshold, float[] respbins, float[] pbins, float[] thetabins, float[] phibins, float effLow, Boolean print){

      float TP=0,FP=0,FN=0;

      int nrespBins=0;
      for(float bin=respbins[1];bin<(respbins[2]-respbins[0]);bin+=respbins[0]){nrespBins++;}
      float[] bincentre_resp=new float[nrespBins];
      float[] TP_resp=new float[nrespBins];
      float[] FP_resp=new float[nrespBins];
      float[] FN_resp=new float[nrespBins];
      for(int i=0;i<nrespBins;i++){
        bincentre_resp[i]=0;
        TP_resp[i]=0;
        FP_resp[i]=0;
        FN_resp[i]=0;
      }

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
  
      Leaf part = new Leaf(32200, 99, "i", 1200);
      Leaf pred_part = new Leaf(32200, 1, "i", 1200);

      /*for priting */
      Leaf trackbank = new Leaf(32000,1,"i",4096);
      Leaf ecalbank = new Leaf(32200,2,"i",4096);
      Leaf htccbank = new Leaf(32200,98,"i",4096);
      Leaf ftofbank = new Leaf(32200,3,"i",4096);
      Bank recpart = r.getBank("REC::Particle");
      Bank reccal = r.getBank("REC::Calorimeter");
      Bank rechtcc = r.getBank("REC::Cherenkov");

      H1F hRespPos = new H1F("e-", 100, 0, 1);
      hRespPos.attr().setLineColor(2);
      hRespPos.attr().setFillColor(2);
      hRespPos.attr().setLineWidth(3);
      hRespPos.attr().setTitleX("Response");
      H1F hRespNeg = new H1F("Not e-", 100, 0, 1);
      hRespNeg.attr().setLineColor(5);
      hRespNeg.attr().setLineWidth(3);
      hRespNeg.attr().setTitleX("Response");
      
      while(r.hasNext() && count<limTotNegEvs){

        r.nextEvent(ev);
        ev.read(part,32200,99);
        ev.read(pred_part,32200,1);

        /* for print */
        if(print){
          ev.read(trackbank,32000,1);
          ev.read(ecalbank,32200,2);
          ev.read(htccbank,32200,98);
          ev.read(ftofbank,32200,3);
          ev.read(recpart);
          ev.read(reccal);
          ev.read(rechtcc);
        }

       
        //same rows in part bank and pred_part bank
        for(int row=0;row<pred_part.getRows();row++){
          short charge = pred_part.getShort(4,row);
          short sector = pred_part.getShort(3,row);
          int recpid=part.getInt(1,row);
          int pindex=part.getShort(0,row);
          float resp=(float)pred_part.getDouble(2,row);
          float px=(float)pred_part.getDouble(6,row);
          float py=(float)pred_part.getDouble(7,row);
          float pz=(float)pred_part.getDouble(8,row);
          double p=Math.sqrt(px*px+py*py+pz*pz);
          double Theta = (Math.acos(pz / p)*(180/Math.PI));// Math.atan2(Math.sqrt(px*px+py*py),pz);
          double Phi = (Math.atan2(py, px)*(180/Math.PI));
          //require neg parts in FD
          //require tracks matched to rec part
          //only estimating el pid metrics not tracking!
          if(sector>0 && charge==-1 && pindex!=999){

            double totHTCC=pred_part.getDouble(30,row)+pred_part.getDouble(31,row)+pred_part.getDouble(32,row);
            if(recpid==11 && totHTCC!=0){ 
              hRespPos.fill(resp);
               /*priting stuff */
               if(print){
                if(resp<0.01){
                  System.out.println("Have one true e- with resp <0.01:");
                  System.out.println("rec part");
                  recpart.show();
                  System.out.println("part ");
                  part.print();
                  System.out.println("pred part ");
                  pred_part.print();
                  System.out.println("htcc ");
                  htccbank.print();
                  System.out.println("track ");
                  trackbank.print();
                  System.out.println("rec cal");
                  reccal.show();
                  System.out.println("rec htcc");
                  rechtcc.show();
                }
              }
            
            }
            else{hRespNeg.fill(resp);}

            int[] mets=checkmets(resp,threshold,recpid,totHTCC);
            TP+=mets[0];
            FP+=mets[1];
            FN+=mets[2];
            int nbin=0;
            for(float bin=respbins[1];bin<(respbins[2]-respbins[0]);bin+=respbins[0]){
              int[] binmets=checkmets(resp,bin,recpid,totHTCC);
              bincentre_resp[nbin]=bin+(respbins[0]/2);
              TP_resp[nbin]+=binmets[0];
              FP_resp[nbin]+=binmets[1];
              FN_resp[nbin]+=binmets[2];
             
              nbin++;
            }

            nbin=0;
            for(float bin=pbins[1];bin<(pbins[2]-pbins[0]);bin+=pbins[0]){
              if(p>=bin && p<(bin+pbins[0])){
                int[] binmets=checkmets(resp,threshold,recpid,totHTCC);
                bincentre_p[nbin]=bin+(pbins[0]/2);
                TP_p[nbin]+=binmets[0];
                FP_p[nbin]+=binmets[1];
                FN_p[nbin]+=binmets[2];
              }
              nbin++;
            }

            nbin=0;
            for(float bin=thetabins[1];bin<(thetabins[2]-thetabins[0]);bin+=thetabins[0]){
              if(Theta>=bin && Theta<(bin+thetabins[0])){
                int[] binmets=checkmets(resp,threshold,recpid,totHTCC);
                bincentre_theta[nbin]=bin+(thetabins[0]/2);
                TP_theta[nbin]+=binmets[0];
                FP_theta[nbin]+=binmets[1];
                FN_theta[nbin]+=binmets[2];
              }
              nbin++;
            }

            nbin=0;
            for(float bin=phibins[1];bin<(phibins[2]-phibins[0]);bin+=phibins[0]){
              if(Phi>=bin && Phi<(bin+phibins[0])){
                int[] binmets=checkmets(resp,threshold,recpid,totHTCC);
                bincentre_phi[nbin]=bin+(phibins[0]/2);
                TP_phi[nbin]+=binmets[0];
                FP_phi[nbin]+=binmets[1];
                FN_phi[nbin]+=binmets[2];
              }
              nbin++;
            }
          }
        }
      }

      TDirectory.export("plots/ElPID0" + ".twig", "/ai/validation/Resp", hRespPos);
      TDirectory.export("plots/ElPID0" + ".twig", "/ai/validation/Resp", hRespNeg);

      System.out.printf("\nTP %f FP %f FN %f \n",TP,FP,FN);

      makeGraphs(nrespBins,effLow,bincentre_resp,TP_resp,FP_resp,FN_resp,"resp","Response","");
      makeGraphs(npBins,effLow,bincentre_p,TP_p,FP_p,FN_p,"P","Momentum","[GeV]");
      makeGraphs(nthetaBins,effLow,bincentre_theta,TP_theta,FP_theta,FN_theta,"theta","Theta","[Degrees]");
      makeGraphs(nphiBins,effLow,bincentre_phi,TP_phi,FP_phi,FN_phi,"phi","Phi","[Degrees]");


      
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
    //read plots in j4shell with eg TwigStudio.browser("plots/ElPID0.twig");
    
    public static void main(String[] args){
        
      System.out.println("----- starting electron PID validator ");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.parse(args);

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
        
      ElPIDValidater dp = new ElPIDValidater();
      dp.process(p.getOption("-in").stringValue(),150000,0.075,respbins,pbins,thetabins,phibins,(float)0.99,false);
  
      
      
        
    }
}
