/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.validation;

import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Bank;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.instarec.utils.DataEntry;
import j4np.instarec.utils.DataList;
import j4np.utils.io.OptionParser;
import twig.data.H1F;
import twig.data.TDirectory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author tyson
 */
public class ElasticScatteringValidator {
    
    public ElasticScatteringValidator(){
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

    public void calcExc(double[] p, double[] el, double[] exc,double beamE, int targPID){
      
      double elE=Math.sqrt(square(el[1])+square(getM(11)));
      double pE=Math.sqrt(square(p[1])+square(getM(targPID)));

      double pTheta=p[2]*(Math.PI/180);
      double elTheta=el[2]*(Math.PI/180);

      double IM = Math.sqrt(square(pE+elE)- ( square(p[6]+el[6]) + square(p[7]+el[7]) + square(p[8]+el[8]) ));
      double px_m = -1.0*(el[6]+p[6]);
      double py_m = -1.0*(el[7]+p[7]);
      double pz_m = beamE - (el[8]+p[8]);
      double p_m=Math.sqrt(square(px_m) + square(py_m)+square(pz_m));
      double pxp_m = px_m/p_m;
      double pyp_m = py_m/p_m;
      double E_m = beamE + getM(targPID) - (elE+pE);
      double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));

      double px_mp = -1.0*(el[6]);
      double py_mp = -1.0*(el[7]);
      double pz_mp = beamE - (el[8]);
      double E_mp = beamE + getM(targPID) - (elE);
      double MM2_p = square(E_mp) - (square(px_mp) + square(py_mp)+square(pz_mp));

      double calcThetaEl=2*Math.atan(Math.cos(pTheta)/(Math.sin(pTheta)*( (beamE/getM(targPID))+1 ) ) );
      double p_i = (getM(targPID)/(1-Math.cos(elTheta)))*(Math.cos(elTheta)+Math.sin(elTheta)*(Math.cos(pTheta)/Math.sin(pTheta))-1);
      double p_f = p_i/(1+(p_i/getM(targPID))*(1-Math.cos(elTheta)));

      exc[0]=IM;
      exc[1]=Math.sqrt(MM2_p);
      exc[2]=Math.sqrt(pxp_m*pxp_m + pyp_m*pyp_m);
      exc[3]=p_i;
      exc[4]=p_f;
      exc[5]=calcThetaEl*(180/Math.PI);
    }

    public void fillRECPart(Bank RECPart, int pindex, double[] part){
      double pz = RECPart.getFloat("pz", pindex);
      double px = RECPart.getFloat("px", pindex);
      double py = RECPart.getFloat("py", pindex);
      double p=Math.sqrt(px*px+py*py+pz*pz);
      double Theta = Math.acos(pz / p)*(180/Math.PI);// Math.atan2(Math.sqrt(px*px+py*py),pz);
      double Phi = Math.atan2(py, px)*(180/Math.PI);
      int pid=RECPart.getInt("pid", pindex);
      int status=RECPart.getInt("status", pindex);
      int charge=RECPart.getInt("charge", pindex);
      part[0]=pid;
      part[1]=p;
      part[2]=Theta;
      part[3]=Phi;
      part[4]=status;
      part[5]=charge;
      part[6]=px;
      part[7]=py;
      part[8]=pz;
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

    public int hasTriggerEl(Bank triggerBank){
      long bits = triggerBank.getLong("trigger", 0);
      int[] L1trigger = convertL1Trigger(bits);
      return L1trigger[0];
    }

    public int hasPredEl(Leaf pred_part, double th){
      int hasEl=0;
      for(int row=0;row<pred_part.getRows();row++){
        short charge = pred_part.getShort(4,row);
        float resp=(float)pred_part.getDouble(2,row);
        if(resp>th && charge==-1){
          hasEl=11;
        }
      }
      return hasEl;
    }

    public int matchTracks(Leaf pred_part, double recpart[],double lim_p,double th){
      int isMatched=-1;
      for(int row=0;row<pred_part.getRows();row++){
        short charge = pred_part.getShort(4,row);
        float px=(float)pred_part.getDouble(6,row);
        float py=(float)pred_part.getDouble(7,row);
        float pz=(float)pred_part.getDouble(8,row);
        float resp=(float)pred_part.getDouble(2,row);
        double p=Math.sqrt(px*px+py*py+pz*pz);
        double Theta = (Math.acos(pz / p)*(180/Math.PI));// Math.atan2(Math.sqrt(px*px+py*py),pz);
        double Phi = (Math.atan2(py, px)*(180/Math.PI));
        double res_onp=Math.abs(recpart[1]-p);
        if(res_onp<(lim_p*recpart[1]) && recpart[5]==charge){ 
          isMatched=211;
          if(resp>th && charge==-1){
            isMatched=11;
          }
        }
      }
      return isMatched;
    }

    public void cleanArr(double[] arr, double length){
      for(int i=0;i<length;i++){arr[i]=0;}
    }

    public int getCalInfo(Bank RECCal, int pindex, double P, double Es[], double Ls[]){
    
      float PCALE=0,ECINE=0,ECOUTE=0;
      int sector=0;
      for (int k = 0; k < RECCal.getRows(); k++) {
        short i = RECCal.getShort("pindex", k);
        int sect=RECCal.getInt("sector", k);
        float lu=RECCal.getFloat("lu",k);
        float lv=RECCal.getFloat("lv",k);
        float lw=RECCal.getFloat("lw",k);
        float energy=RECCal.getFloat("energy",k);
        byte layer=RECCal.getByte("layer",k);
        if (i == pindex ) {
          sector=sect;
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
      return sector;
    }

    public Boolean passFid(double[] Ls){
      if(Ls[0]>9 && Ls[1]>9 && Ls[2]>9){
        return true;
      } else{
        return false;
      }
    }

    public Double getHTCCSectorSum(Leaf htcc,int sector){
      double sum=0;
      if(sector!=0){
        for(int i=1;i<9;i++){
          sum+=htcc.getDouble(i, sector-1);
          //System.out.printf("i %d sector %d htcc sector %d htcc val %f sum %f\n",i,sector,htcc.getInt(0, sector-1),htcc.getDouble(i, sector-1),sum);
        }
      }
      return sum;
    }


    public void process(String file, int limEvs,String endName, double threshold, int targPID, double beamE, Boolean requireMatchedEl, Boolean print, int desired_sector, int... tags){
      double lim_p_res=0.2; //percentage

      //output training sample size per sector
      int count = 0;

      HipoReader r = new HipoReader();
      Event ev = new Event();
      r.setTags(11);
      r.open(file);
      

      /*for priting */
      Bank recpart = r.getBank("REC::Particle");
      Bank reccal = r.getBank("REC::Calorimeter");
      Bank triggerbank = r.getBank("RUN::config");

      Leaf pred_part = new Leaf(32200, 1, "i", 1200);
      Leaf htcc = new Leaf(32200,98, "i", 4096);

      int nVarsPart=9,nVarsExc=6;
      // pid, P, Theta, Phi
      double[] p=new double[nVarsPart];
      double[] el=new double[nVarsPart];
      double[] exc=new double[nVarsExc];
      double[] elEs=new double[7];
      double[] elLs=new double[9];

      double targMass=getM(targPID);

      H1F hMM_nop = new H1F("Missing Recoil Mass", 100,0.5,3.5);
      hMM_nop.attr().setLineColor(2);
      hMM_nop.attr().setLineWidth(3);
      hMM_nop.attr().setTitleX("Mass [GeV]");

      H1F hdPhi = new H1F("dPhi", 100,165,195);
      hdPhi.attr().setLineColor(2);
      hdPhi.attr().setLineWidth(3);
      hdPhi.attr().setTitleX("dPhi [Degree]");

      H1F hPtP = new H1F("PtP", 100,0.,0.1);
      hPtP.attr().setLineColor(2);
      hPtP.attr().setLineWidth(3);
      hPtP.attr().setTitleX("PtP [GeV]");

      H1F hdThetaEl = new H1F("dThetaEl", 100,-20,20);
      hdThetaEl.attr().setLineColor(2);
      hdThetaEl.attr().setLineWidth(3);
      hdThetaEl.attr().setTitleX("dThetaEl [Degree]");

      H1F hdBeamE = new H1F("dBeamE", 100,-5,5);
      hdBeamE.attr().setLineColor(2);
      hdBeamE.attr().setLineWidth(3);
      hdBeamE.attr().setTitleX("dBeamE [GeV]");

      H1F hdP = new H1F("With L1 e-", 100,-5,5);
      hdP.attr().setLineColor(2);
      hdP.attr().setLineWidth(3);
      hdP.attr().setTitleX("dP [GeV]");

      H1F hdP_wInstaEl = new H1F("With L1 e- & Online e-", 100,-5,5);
      hdP_wInstaEl.attr().setLineColor(5);
      hdP_wInstaEl.attr().setLineWidth(3);
      hdP_wInstaEl.attr().setTitleX("dP [GeV]");

      H1F hdP_wOffEl = new H1F("With L1 e- & Offline e-", 100,-5,5);
      hdP_wOffEl.attr().setLineColor(1);
      hdP_wOffEl.attr().setLineWidth(3);
      hdP_wOffEl.attr().setTitleX("dP [GeV]");

      H1F hdP_radel = new H1F("With L1 e- (Radiative Elastic)", 100,-5,5);
      hdP_radel.attr().setLineColor(2);
      hdP_radel.attr().setLineWidth(3);
      hdP_radel.attr().setTitleX("dP [GeV]");

      H1F hdP_wInstaEl_radel = new H1F("With L1 e- & Online e- (Radiative Elastic)", 100,-5,5);
      hdP_wInstaEl_radel.attr().setLineColor(5);
      hdP_wInstaEl_radel.attr().setLineWidth(3);
      hdP_wInstaEl_radel.attr().setTitleX("dP [GeV]");

      H1F hdP_wOffEl_radel = new H1F("With L1 e- & Offline e- (Radiative Elastic)", 100,-5,5);
      hdP_wOffEl_radel.attr().setLineColor(1);
      hdP_wOffEl_radel.attr().setLineWidth(3);
      hdP_wOffEl_radel.attr().setTitleX("dP [GeV]");

      H1F hP = new H1F("Elastic Scattering", 100,0,11);
      hP.attr().setLineColor(2);
      hP.attr().setLineWidth(3);
      hP.attr().setTitleX("P [GeV]");

      H1F hP_radel = new H1F("Radiative Elastic Scattering", 100,0,11);
      hP_radel.attr().setLineColor(6);
      hP_radel.attr().setLineWidth(3);
      hP_radel.attr().setTitleX("P [GeV]");

      //don't count if set limEvs=-1
      if(limEvs==-1){count=-2;}
      
      while(r.hasNext() && count<limEvs){

        r.nextEvent(ev);
        ev.read(pred_part,32200,1);
        ev.read(htcc,32200,98);
        ev.read(recpart);
        ev.read(reccal);
        ev.read(triggerbank);

        for(int row=0;row<recpart.getRows();row++){
          for (int row2 = 0; row2 < recpart.getRows(); row2++) {

            int desired_sector_e=desired_sector;

            cleanArr(p,nVarsPart);
            cleanArr(el,nVarsPart);
            cleanArr(exc,nVarsExc);
            cleanArr(elEs,7);
            cleanArr(elLs,9);

            fillRECPart(recpart,row2,p);
            fillRECPart(recpart,row,el);
            calcExc(p, el, exc,beamE,targPID);

            int hasTrig=hasTriggerEl(triggerbank);
            int hasPredEl=hasPredEl(pred_part, threshold);

            int elmatch=-1;
            if(requireMatchedEl){
              elmatch=matchTracks(pred_part,el,lim_p_res,threshold);
            } else{
              elmatch=hasPredEl;
            }



            //System.out.println("\n New comb");
            //System.out.printf("p1 pid %f p %f theta %f phi %f status %f charge %f\n",p[0],p[1],p[2],p[3],p[4],p[5]);
            //System.out.printf("p2 pid %f p %f theta %f phi %f status %f charge %f\n",el[0],el[1],el[2],el[3],el[4],el[5]);


            int sector=getCalInfo(reccal, row, el[1], elEs, elLs);
            Boolean goodEl=false;
            double totHTCC=getHTCCSectorSum(htcc,sector); //&& totHTCC!=0

            //System.out.printf("sector %d totHTCC %f \n",sector,totHTCC);
            //htcc.print();

            if(desired_sector==0){
              desired_sector_e=sector;
            }

            if(elEs[0]>0.06 && elEs[3]>0.15  && totHTCC!=0 && passFid(elLs) && desired_sector_e==sector){ //&& elEs[3]>0.15 && passFid(elLs)
              goodEl=true;
            }

            if( Math.abs(p[4])>4000 && p[5]==1 && p[0]==targPID){
              if(Math.abs(el[4])>=2000 && Math.abs(el[4])<4000 && el[5]==-1 && elmatch!=-1 && goodEl){

                //System.out.println("\n past cuts");
                //System.out.printf("p1 pid %f p %f theta %f phi %f status %f charge %f\n",p[0],p[1],p[2],p[3],p[4],p[5]);
                //System.out.printf("p2 pid %f p %f theta %f phi %f status %f charge %f\n",el[0],el[1],el[2],el[3],el[4],el[5]);
                if(hasTrig!=0){
                  double phiDif=p[3]-el[3];
                  hdPhi.fill(phiDif);
                  if(phiDif>175 && phiDif<185 ){
                    //don't count if set limEvs=-1
                    if(limEvs!=-1){
                      count++;
                    }
                    hMM_nop.fill(exc[1]);
                    hPtP.fill(exc[2]);

                    //elastic
                    if(exc[1]>(targMass-0.6) && exc[1]<(targMass+0.6)){
                      hdBeamE.fill(exc[3]-beamE);
                      hdP.fill(exc[4]-el[1]);
                      hP.fill(el[1]);
                      hdThetaEl.fill(exc[5]-el[2]);
                      if(elmatch==11){
                        hdP_wInstaEl.fill(exc[4]-el[1]);
                      }
                      if(el[0]==11){
                        hdP_wOffEl.fill(exc[4]-el[1]);
                      }
                    }
                    //radiative elastic
                    if(exc[2]<0.025){
                      hdP_radel.fill(exc[4]-el[1]);
                      hP_radel.fill(el[1]);
                      if(elmatch==11){ //elmatch, hasPredEl
                        hdP_wInstaEl_radel.fill(exc[4]-el[1]);
                      }
                      if(el[0]==11){
                        hdP_wOffEl_radel.fill(exc[4]-el[1]);
                      }

                    }

                  }//good phi dif

                }//has e- trigger bit
              }//good el candidate
            }//good recoil nucleon candidate

          }//loop over nucleon candidates
        }//loop over electron candidates

      }//read file

      if(desired_sector==0){
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdPhi);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hMM_nop);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hPtP);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdBeamE);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdThetaEl);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_wInstaEl);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_wOffEl);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_radel);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_wInstaEl_radel);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_wOffEl_radel);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hP);
        TDirectory.export("plots/elasticScattering"+endName+".twig","/ai/validation/tPID"+String.valueOf(targPID),hP_radel);
      } else {
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdPhi);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hMM_nop);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hPtP);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdBeamE);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdThetaEl);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_wInstaEl);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_wOffEl);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_radel);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_wInstaEl_radel);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hdP_wOffEl_radel);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hP);
        TDirectory.export("plots/elasticScattering"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation/tPID"+String.valueOf(targPID),hP_radel);
      }
      
    }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar -in w.h5 
    //read plots in j4shell with eg TwigStudio.browser("plots/Tracking.twig");
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting elastic scattering validator ");
      OptionParser p = new OptionParser();
      p.addRequired("-in", "input name");
      p.parse(args);

      String endName="_wCuts_th0p05_phiCut5_fromTagger"; //_phiCut5 eg 175-185, 10 otherwise

      double resp_threshold=0.05;
      double beamE=10.5473;
        
      ElasticScatteringValidator dp = new ElasticScatteringValidator();
      dp.process(p.getOption("-in").stringValue(),-1,endName+"_onProton",resp_threshold,2212,beamE,true,false,0,11);
      //dp.process(p.getOption("-in").stringValue(),-1,endName+"_onDeuteron",resp_threshold,45,beamE,false,false,0);
        
    }

}