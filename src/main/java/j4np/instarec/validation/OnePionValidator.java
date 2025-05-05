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
import twig.data.BarChartBuilder;
import twig.data.DataGroup;
import twig.data.H1F;
import twig.data.H2F;
import twig.data.TDirectory;
import twig.graphics.TGCanvas;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

/**
 *
 * @author tyson
 */
public class OnePionValidator {
    
    public OnePionValidator(){
    }

    public void fillHisto(Map<String, Integer> names,Vector<H1F> histos,String name,double value){
      histos.get(names.get(name)).fill(value);
    }

    public void fillHisto2D(Map<String, Integer> names,Vector<H2F> histos,String name,double value,double value2){
      histos.get(names.get(name)).fill(value,value2);
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

    public void calcExc(double[] pip, double[] el, double[] exc,double beamE){
      
      double elE=Math.sqrt(square(el[1])+square(getM(11)));
      double pipE=Math.sqrt(square(pip[1])+square(getM(211)));
      double pM=getM(2212);

      double IM = Math.sqrt(square(pipE+elE)- ( square(pip[6]+el[6]) + square(pip[7]+el[7]) + square(pip[8]+el[8]) ));
      double px_m = -1.0*(el[6]+pip[6]);
      double py_m = -1.0*(el[7]+pip[7]);
      double pz_m = beamE - (el[8]+pip[8]);
      double p_m=Math.sqrt(square(px_m) + square(py_m)+square(pz_m));
      double pxp_m = px_m/p_m;
      double pyp_m = py_m/p_m;
      double E_m = beamE + getM(2212) - (elE+pipE);
      double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));

      exc[0]=IM;
      exc[1]=Math.sqrt(MM2);
      exc[2]=Math.sqrt(pxp_m*pxp_m + pyp_m*pyp_m);
      exc[3]=(MM2-pM*pM)/(2*pM);
    }

    public void fillRECPart(Bank RECPart, Bank trackBank, int pindex, double[] part){
      double pz = RECPart.getFloat("pz", pindex);
      double vz = RECPart.getFloat("vz", pindex);
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
      part[10]=0;
      part[13]=vz;
      //part 9 is track sector, 11 is track chi^2, 12 is nb of track superlayers
      fillTrackInfo(pindex,trackBank,part);
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
      int bestIsMatched=-1;
      double bestRes_onp=999;
      double sumHTCC=0;
      for(int row=0;row<pred_part.getRows();row++){
        short charge = pred_part.getShort(4,row);
        short sector = pred_part.getShort(3,row);
        float px=(float)pred_part.getDouble(6,row);
        float py=(float)pred_part.getDouble(7,row);
        float pz=(float)pred_part.getDouble(8,row);
        float resp=(float)pred_part.getDouble(2,row);
        int pid=pred_part.getInt(1,row);
        double p=Math.sqrt(px*px+py*py+pz*pz);
        double Theta = (Math.acos(pz / p)*(180/Math.PI));// Math.atan2(Math.sqrt(px*px+py*py),pz);
        double Phi = (Math.atan2(py, px)*(180/Math.PI));
        double res_onp=Math.abs(recpart[1]-p);
        if(res_onp<Math.abs(lim_p*recpart[1]) && recpart[5]==charge && recpart[9]==sector){ 
          isMatched=charge*211;
          if(resp>th && charge==-1){ //resp>th
            isMatched=11;
          }
          if(res_onp<bestRes_onp){
            bestRes_onp=res_onp;
            bestIsMatched=isMatched;
            sumHTCC=pred_part.getDouble(30,row)+pred_part.getDouble(31,row)+pred_part.getDouble(32,row);
          }
        }
      }
      recpart[10]=sumHTCC;
      return bestIsMatched;
    }

    public void cleanArr(double[] arr, double length){
      for(int i=0;i<length;i++){arr[i]=0;}
    }

    public void copyArr(double[] arr,double[] arr_cp, double length){
      for(int i=0;i<length;i++){arr_cp[i]=arr[i];}
    }

    public int getCalInfo(Bank RECCal, int pindex, double P, double Es[], double Ls[]){
    
      float PCALE=0,ECINE=0,ECOUTE=0,X=0,Y=0;
      int sector=0;
      for (int k = 0; k < RECCal.getRows(); k++) {
        short i = RECCal.getShort("pindex", k);
        int sect=RECCal.getInt("sector", k);
        float lu=RECCal.getFloat("lu",k);
        float lv=RECCal.getFloat("lv",k);
        float lw=RECCal.getFloat("lw",k);
        float x=RECCal.getFloat("x",k);
        float y=RECCal.getFloat("y",k);
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
            X=x;
            Y=y;
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
      Es[7]=X;
      Es[8]=Y;
      return sector;
    }

    public void fillTrackInfo(int pindext, Bank trackBank,double[] part){
      int s=0,nSL=6;
      float chi2=9999;
      for(int i=0;i<trackBank.getRows();i++){
        int sector=trackBank.getInt("sector",i);
        int pindex=trackBank.getInt("pindex",i);
        float c2=trackBank.getFloat("chi2", i);
        if(pindex==pindext){
          s=sector;
          chi2=c2;
          if( (trackBank.getInt("status", i) & 0b101010101010)!=0){
            nSL=5;
          }
        }
      }
      part[9]=s;
      part[11]=chi2;
      part[12]=nSL;
    }

    public Boolean passFid(double[] Ls, double part[],int str){
      //is track chi^2 && part11]<350 && Math.abs(part13])<20 && part12]==6 && oneOfEcalHTCC==1
      //RGA vz && part[13]<20 && part[13]>(-13)
      //tight: >14cm ortherwise 9cm

      boolean vz=false,looseCal=false,tightCal=false;

      //for str 0:
      if(part[13]<20 && part[13]>(-13)){
        vz=true;
      }

      if(Ls[0]>9 && Ls[1]>9 && Ls[2]>9){
        looseCal=true;
        if(Ls[0]>14 && Ls[1]>14 && Ls[2]>14){
          tightCal=true;
        }
      }

      if(str==0){
        if(vz){
          return true;
        } else{
          return false;
        }
      } else if (str==1){
        if(vz && looseCal){
          return true;
        } else{
          return false;
        }
      }  else if (str==2){
        if(vz && tightCal){
          return true;
        } else{
          return false;
        }
      } else {
        return true;
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


    public void process(String file, int limEvs,String endName, double threshold,double beamE, int reqFids, int desired_sector){
      
      double lim_p_res=0.2; //percentage

      //output training sample size per sector
      int count = 0, nEl=0, nNotEl=0, nAllInstaEl=0;

      HipoReader r = new HipoReader();
      Event ev = new Event();
      r.open(file);
      

      /*for printing */
      Bank recpart = r.getBank("REC::Particle");
      Bank reccal = r.getBank("REC::Calorimeter");
      Bank rectrack = r.getBank("REC::Track");
      Bank triggerbank = r.getBank("RUN::config");

      Leaf pred_part = new Leaf(32, 3, "i", 1200);

      int nVarsPart=14,nVarsExc=4;
      double[] pip=new double[nVarsPart];
      double[] el=new double[nVarsPart];
      double[] exc=new double[nVarsExc];
      double[] elEs=new double[9];
      double[] elLs=new double[9];
      double[] pipEs=new double[9];
      double[] pipLs=new double[9];

      double targMass=getM(2212);

      Vector<H1F> histos = new Vector<>();
      Map<String, Integer> names = new HashMap<>();

      Vector<H2F> histos2D = new Vector<>();
      Map<String, Integer> names2D = new HashMap<>();

      int nHistos=0,nHistos2D=0;
      for(int i=0;i<3;i++){
        String instapidelstring="";
        String instapidelfilestring="";
        if(i==1){
          instapidelstring="!";
          instapidelfilestring="_noInstaEl";
        } else if(i==2){
          instapidelstring="with & without";
          instapidelfilestring="_all";
        }
        
        H1F hNCal = new H1F("Calo: Offline e- PID & "+instapidelstring+" Online e- PID", 6,-1.5,4.5);
        hNCal.attr().setLineColor(5);
        hNCal.attr().setLineWidth(3);
        hNCal.attr().setTitleX("Nb of Calorimeter Layers Hit");
        names.put("hNCal"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hNCal);

        H1F hNCal_w = new H1F("Calo: Offline ! e- PID & "+instapidelstring+" Online e- PID", 6,-1.5,4.5);
        hNCal_w.attr().setLineColor(5);
        hNCal_w.attr().setLineWidth(3);
        hNCal_w.attr().setTitleX("Nb of Calorimeter Layers Hit");
        names.put("hNCal_w"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hNCal_w);

        H1F hNCal_allInstaEl = new H1F("Calo: "+instapidelstring+" Online e- PID", 6,-1.5,4.5);
        hNCal_allInstaEl.attr().setLineColor(5);
        hNCal_allInstaEl.attr().setLineWidth(3);
        hNCal_allInstaEl.attr().setTitleX("Nb of Calorimeter Layers Hit");
        names.put("hNCal_allInstaEl"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hNCal_allInstaEl);

        H1F hPCAL_allInstaEl = new H1F("PCAL: "+instapidelstring+" Online e- PID", 4,-1.5,2.5);
        hPCAL_allInstaEl.attr().setLineColor(5);
        hPCAL_allInstaEl.attr().setLineWidth(3);
        hPCAL_allInstaEl.attr().setTitleX("PCAL is Hit");
        names.put("hPCAL_allInstaEl"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hPCAL_allInstaEl);

        H1F hECIN_allInstaEl = new H1F("ECIN: "+instapidelstring+" Online e- PID", 4,-1.5,2.5);
        hECIN_allInstaEl.attr().setLineColor(5);
        hECIN_allInstaEl.attr().setLineWidth(3);
        hECIN_allInstaEl.attr().setTitleX("ECIN is Hit");
        names.put("hECIN_allInstaEl"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hECIN_allInstaEl);

        H1F hECOUT_allInstaEl = new H1F("ECOUT: "+instapidelstring+" Online e- PID", 4,-1.5,2.5);
        hECOUT_allInstaEl.attr().setLineColor(5);
        hECOUT_allInstaEl.attr().setLineWidth(3);
        hECOUT_allInstaEl.attr().setTitleX("ECOUT is Hit");
        names.put("hECOUT_allInstaEl"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hECOUT_allInstaEl);

        H1F hHTCC_allInstaEl = new H1F("HTCC: "+instapidelstring+" Online e- PID", 4,-1.5,2.5);
        hHTCC_allInstaEl.attr().setLineColor(5);
        hHTCC_allInstaEl.attr().setLineWidth(3);
        hHTCC_allInstaEl.attr().setTitleX("HTCC is Hit");
        names.put("hHTCC_allInstaEl"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hHTCC_allInstaEl);

        H1F hMM = new H1F("Mx: Offline e- PID & "+instapidelstring+" Online e- PID", 100,0,3.5);
        hMM.attr().setLineColor(5);
        hMM.attr().setLineWidth(3);
        hMM.attr().setTitleX("Mx(e'#pi^+) [GeV]");
        names.put("hMM"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hMM);

        H1F hMM_w = new H1F("Mx: Offline ! e- PID & "+instapidelstring+" Online e- PID", 100,0,3.5);
        hMM_w.attr().setLineColor(1);
        hMM_w.attr().setLineWidth(3);
        hMM_w.attr().setTitleX("Mx(e'#pi^+) [GeV]");
        names.put("hMM_w"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hMM_w);

        H1F hMM_allInstaEl = new H1F("Mx: "+instapidelstring+" Online e- PID", 100,0,3.5);
        hMM_allInstaEl.attr().setLineColor(2);
        hMM_allInstaEl.attr().setLineWidth(3);
        hMM_allInstaEl.attr().setTitleX("Mx(e'#pi^+) [GeV]");
        names.put("hMM_allInstaEl"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hMM_allInstaEl);

        H1F hIM = new H1F("M: Offline e- PID & "+instapidelstring+" Online e- PID", 100,0,3.5);
        hIM.attr().setLineColor(5);
        hIM.attr().setLineWidth(3);
        hIM.attr().setTitleX("M(#pi^+e^-) [GeV]");
        names.put("hIM"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hIM);

        H1F hIM_w = new H1F("M: Offline ! e- PID & "+instapidelstring+" Online e- PID", 100,0,3.5);
        hIM_w.attr().setLineColor(1);
        hIM_w.attr().setLineWidth(3);
        hIM_w.attr().setTitleX("M(#pi^+e^-) [GeV]");
        names.put("hIM_w"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hIM_w);

        H1F hIM_allInstaEl = new H1F("M: "+instapidelstring+" Online e- PID", 100,0,3.5);
        hIM_allInstaEl.attr().setLineColor(2);
        hIM_allInstaEl.attr().setLineWidth(3);
        hIM_allInstaEl.attr().setTitleX("M(#pi^+e^-) [GeV]");
        names.put("hIM_allInstaEl"+instapidelstring,nHistos);
        nHistos++;
        histos.add(hIM_allInstaEl);

        H2F helPTheta = new H2F("e^- #theta P: Offline e- PID & "+instapidelstring+" Online e- PID",50,0,10,50,0,40);
        helPTheta.attr().setTitleX("P [GeV]");
        helPTheta.attr().setTitleY("#theta [Deg]");
        names2D.put("helPTheta"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helPTheta);

        H2F helPTheta_w = new H2F("e^- #theta P: Offline ! e- PID & "+instapidelstring+" Online e- PID",50,0,10,50,0,40);
        helPTheta_w.attr().setTitleX("P [GeV]");
        helPTheta_w.attr().setTitleY("#theta [Deg]");
        names2D.put("helPTheta_w"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helPTheta_w);

        H2F hpipPTheta_allInstaEl = new H2F("#pi^+ #theta P:  "+instapidelstring+" Online e- PID ",50,0,10,50,0,40);
        hpipPTheta_allInstaEl.attr().setTitleX("P [GeV]");
        hpipPTheta_allInstaEl.attr().setTitleY("#theta [Deg]");
        names2D.put("hpipPTheta_allInstaEl"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(hpipPTheta_allInstaEl);


        H2F helPTheta_allInstaEl = new H2F("e^- #theta P:  "+instapidelstring+" Online e- PID",50,0,10,50,0,40);
        helPTheta_allInstaEl.attr().setTitleX("P [GeV]");
        helPTheta_allInstaEl.attr().setTitleY("#theta [Deg]");
        names2D.put("helPTheta_allInstaEl"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helPTheta_allInstaEl);


        H2F helXY = new H2F("e^- Y vs X: Offline e- PID & "+instapidelstring+" Online e- PID",100,-250,250,100,-250,250);
        helXY.attr().setTitleX("X [cm]");
        helXY.attr().setTitleY("Y [cm]");
        names2D.put("helXY"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helXY);

        H2F helXY_w = new H2F("e^- Y vs X: Offline ! e- PID & "+instapidelstring+" Online e- PID",100,-250,250,100,-250,250);
        helXY_w.attr().setTitleX("X [cm]");
        helXY_w.attr().setTitleY("Y [cm]");
        names2D.put("helXY_w"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helXY_w);

        H2F helXY_allInstaEl = new H2F("e^- Y vs X:  "+instapidelstring+" Online e- PID",100,-250,250,100,-250,250);
        helXY_allInstaEl.attr().setTitleX("X [cm]");
        helXY_allInstaEl.attr().setTitleY("Y [cm]");
        names2D.put("helXY_allInstaEl"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helXY_allInstaEl);

        H2F helPSF = new H2F("e^- SF vs P: Offline e- PID & "+instapidelstring+" Online e- PID",50,0,10,50,0,0.5);
        helPSF.attr().setTitleX("P [GeV]");
        helPSF.attr().setTitleY("Sampling Fraction");
        names2D.put("helPSF"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helPSF);

        H2F helPSF_w = new H2F("e^- SF vs P: Offline ! e- PID & "+instapidelstring+" Online e- PID",50,0,10,50,0,0.5);
        helPSF_w.attr().setTitleX("P [GeV]");
        helPSF_w.attr().setTitleY("Sampling Fraction");
        names2D.put("helPSF_w"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helPSF_w);

        H2F helPSF_allInstaEl = new H2F("e^- SF vs P:  "+instapidelstring+" Online e- PID",50,0,10,50,0,0.5);
        helPSF_allInstaEl.attr().setTitleX("P [GeV]");
        helPSF_allInstaEl.attr().setTitleY("Sampling Fraction");
        names2D.put("helPSF_allInstaEl"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helPSF_allInstaEl);

        H2F helLVSF = new H2F("e^- SF vs LV: Offline e- PID & "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
        helLVSF.attr().setTitleX("LV [cm]");
        helLVSF.attr().setTitleY("Sampling Fraction");
        names2D.put("helLVSF"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helLVSF);

        H2F helLVSF_w = new H2F("e^- SF vs LV: Offline ! e- PID & "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
        helLVSF_w.attr().setTitleX("LV [cm]");
        helLVSF_w.attr().setTitleY("Sampling Fraction");
        names2D.put("helLVSF_w"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helLVSF_w);

        H2F helLVSF_allInstaEl = new H2F("e^- SF vs LV:  "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
        helLVSF_allInstaEl.attr().setTitleX("LV [cm]");
        helLVSF_allInstaEl.attr().setTitleY("Sampling Fraction");
        names2D.put("helLVSF_allInstaEl"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helLVSF_allInstaEl);

        H2F helLWSF = new H2F("e^- SF vs LW: Offline e- PID & "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
        helLWSF.attr().setTitleX("LW [cm]");
        helLWSF.attr().setTitleY("Sampling Fraction");
        names2D.put("helLWSF"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helLWSF);

        H2F helLWSF_w = new H2F("e^- SF vs LW: Offline ! e- PID & "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
        helLWSF_w.attr().setTitleX("LW [cm]");
        helLWSF_w.attr().setTitleY("Sampling Fraction");
        names2D.put("helLWSF_w"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helLWSF_w);

        H2F helLWSF_allInstaEl = new H2F("e^- SF vs LW:  "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
        helLWSF_allInstaEl.attr().setTitleX("LW [cm]");
        helLWSF_allInstaEl.attr().setTitleY("Sampling Fraction");
        names2D.put("helLWSF_allInstaEl"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helLWSF_allInstaEl);

      }

      

      //don't count if set limEvs=-1
      if(limEvs==-1){count=-2;}
      
      while(r.hasNext() && count<limEvs){
        //don't count if set limEvs=-1
        if(limEvs!=-1){
          count++;
        }

        r.nextEvent(ev);
        //Leaf pred_part = ev.readLeaf(1,12,32,3);
        ev.read(pred_part,32,3);
        ev.read(recpart);
        ev.read(reccal);
        ev.read(rectrack);
        ev.read(triggerbank);

        int hasTrig=hasTriggerEl(triggerbank);
        

        //System.out.println("\n\nNew Event");
        //pred_part.print();
        //reccal.show();

        if(pred_part.getSize()>0){

          for(int row=0;row<recpart.getRows();row++){
            for(int row3=0;row3<recpart.getRows();row3++){

              cleanArr(pip,nVarsPart);
              cleanArr(el,nVarsPart);
              cleanArr(exc,nVarsExc);
              cleanArr(elEs,9);
              cleanArr(elLs,9);
              cleanArr(pipEs,9);
              cleanArr(pipLs,9);

              fillRECPart(recpart,rectrack,row,el);
              getCalInfo(reccal, row, el[1], elEs, elLs);
              int matchel=matchTracks(pred_part,el,lim_p_res,threshold);
              fillRECPart(recpart,rectrack,row3,pip);
              getCalInfo(reccal, row3, pip[1], pipEs, pipLs);
              int matchpip=matchTracks(pred_part,pip,lim_p_res,threshold);
                

              //System.out.printf("p pid %f px %f py %f pz %f status %f sector %f charge %f match %d\n",part[0],part[6],part[7],part[8],part[4],part[9],part[5],match);

              int desired_sector_e=desired_sector;
              int sector=(int)el[9];
              //use all sectors if desired sector is 0
              if(desired_sector==0){
                desired_sector_e=sector;
              }

              Boolean fid=passFid(elLs,el,reqFids);

              int hasEl=0, hasRECEl=0, hasElCandi=0,hasPip=0;
              int elNCal=0, hasPCAL=0,hasECIN=0,hasECOUT=0,hasHTCC=0,oneOfEcalHTCC=0,ecalHTCC=0;
              if(elEs[0]>0.01){elNCal++;hasPCAL=1;}
              if(elEs[1]>0.01){elNCal++;hasECIN=1;}
              if(elEs[2]>0.01){elNCal++;hasECOUT=1;}
              if(el[10]>0.0){hasHTCC=1;}
              if(hasHTCC>0 || elNCal>0){oneOfEcalHTCC=1;}
              if(hasHTCC>0 && elNCal>0){ecalHTCC=1;}


              if(pip[5]==1 && matchpip!=-1 && pip[0]==211){
                hasPip=1;
              }

              //&& matchel!=-1 el[11] 
              if(el[5]==-1 && desired_sector_e==sector && fid && matchel!=-1 && ecalHTCC==1){ 
                hasElCandi=1;
              }

              if(hasElCandi==1 && hasPip==1 && hasTrig==1){

                calcExc(pip, el, exc,beamE);
                if(el[1]>2 && exc[1]>0.2){ //exc[1]>0.2 && el[1]>2

                  //System.out.println("\n\nFound particles");
                  //System.out.printf("el pid %f px %f py %f pz %f status %f sector %f charge %f\n",el[0],el[6],el[7],el[8],el[4],el[9],el[5]);
                  //System.out.printf("pi- pid %f px %f py %f pz %f status %f sector %f charge %f\n",pim[0],pim[6],pim[7],pim[8],pim[4],pim[9],pim[5]);
                  //System.out.printf("pi+ pid %f px %f py %f pz %f status %f sector %f charge %f\n",pip[0],pip[6],pip[7],pip[8],pip[4],pip[9],pip[5]);

                  for(int j=0;j<2;j++){
                    String instapidelstring="with & without";
                    if(j==1){
                      instapidelstring="!";
                      if(matchel==11){
                        instapidelstring="";
                      }
                    }
                    fillHisto(names,histos,"hIM_allInstaEl"+instapidelstring,exc[0]);
                    fillHisto(names,histos,"hMM_allInstaEl"+instapidelstring,exc[1]);


                    fillHisto2D(names2D,histos2D,"helPTheta_allInstaEl"+instapidelstring,el[1],el[2]);
                    fillHisto2D(names2D,histos2D,"helPSF_allInstaEl"+instapidelstring,el[1],elEs[3]);
                    fillHisto2D(names2D,histos2D,"helXY_allInstaEl"+instapidelstring,elEs[7],elEs[8]);
                    fillHisto2D(names2D,histos2D,"hpipPTheta_allInstaEl"+instapidelstring,pip[1],pip[2]);
                    fillHisto2D(names2D,histos2D,"helLVSF_allInstaEl"+instapidelstring,elLs[1],elEs[3]);
                    fillHisto2D(names2D,histos2D,"helLWSF_allInstaEl"+instapidelstring,elLs[2],elEs[3]);

                    fillHisto(names,histos,"hNCal_allInstaEl"+instapidelstring,elNCal);
                    fillHisto(names,histos,"hPCAL_allInstaEl"+instapidelstring,hasPCAL);
                    fillHisto(names,histos,"hECIN_allInstaEl"+instapidelstring,hasECIN);
                    fillHisto(names,histos,"hECOUT_allInstaEl"+instapidelstring,hasECOUT);
                    fillHisto(names,histos,"hHTCC_allInstaEl"+instapidelstring,hasHTCC);
                      
                    if(el[0]==11){
                      fillHisto(names,histos,"hIM"+instapidelstring,exc[0]);
                      fillHisto(names,histos,"hMM"+instapidelstring,exc[1]);
                      fillHisto2D(names2D,histos2D,"helPTheta"+instapidelstring,el[1],el[2]);
                      fillHisto2D(names2D,histos2D,"helPSF"+instapidelstring,el[1],elEs[3]);
                      fillHisto2D(names2D,histos2D,"helXY"+instapidelstring,elEs[7],elEs[8]);
                      fillHisto2D(names2D,histos2D,"helLVSF"+instapidelstring,elLs[1],elEs[3]);
                      fillHisto2D(names2D,histos2D,"helLWSF"+instapidelstring,elLs[2],elEs[3]);

                      fillHisto(names,histos,"hNCal"+instapidelstring,elNCal);
                    } else {
                      fillHisto(names,histos,"hIM_w"+instapidelstring,exc[0]);
                      fillHisto(names,histos,"hMM_w"+instapidelstring,exc[1]);
                      fillHisto2D(names2D,histos2D,"helPTheta_w"+instapidelstring,el[1],el[2]);
                      fillHisto2D(names2D,histos2D,"helPSF_w"+instapidelstring,el[1],elEs[3]);
                      fillHisto2D(names2D,histos2D,"helXY_w"+instapidelstring,elEs[7],elEs[8]);
                      fillHisto2D(names2D,histos2D,"helLVSF_w"+instapidelstring,elLs[1],elEs[3]);
                      fillHisto2D(names2D,histos2D,"helLWSF_w"+instapidelstring,elLs[2],elEs[3]);

                      fillHisto(names,histos,"hNCal_w"+instapidelstring,elNCal);
                    }  
                  }
                }//fill histos for both with and without online el first then fill with or without
              } //right number of candidates
            }//loop over rec parts 3
          }//loop over rec parts 1
        }//check bank not empty

      }//read file
  
      for (Map.Entry<String, Integer> entry : names.entrySet()) {
        String name=entry.getKey();
        int index=entry.getValue();

        String instapidelfilestring="";
        if(name.contains("!")){instapidelfilestring="_noInstaEl";}
        if(name.contains("with & without")){instapidelfilestring="_all";}

        if(desired_sector==0){
          TDirectory.export("plots_rga/OnePion"+endName+".twig","/ai/validation"+instapidelfilestring,histos.get(index));
        } else {
          TDirectory.export("plots_rga/OnePion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,histos.get(index));
        }
      }

      for (Map.Entry<String, Integer> entry : names2D.entrySet()) {
        String name=entry.getKey();
        int index=entry.getValue();

        String instapidelfilestring="";
        if(name.contains("!")){instapidelfilestring="_noInstaEl";}
        if(name.contains("with & without")){instapidelfilestring="_all";}

        if(desired_sector==0){
          TDirectory.export("plots_rga/OnePion"+endName+".twig","/ai/validation"+instapidelfilestring,histos2D.get(index));
        } else {
          TDirectory.export("plots_rga/OnePion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,histos2D.get(index));
        }
      }
      
    }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar
    //read plots in j4shell with eg TwigStudio.browser("plots/OnePion.twig")
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting one pion validator ");
      //String fName_11="/Users/tyson/data_repo/trigger_data/rga/instarec_skims/cooked_11.h5";
      //String fName_211="/Users/tyson/data_repo/trigger_data/rga/instarec_skims/cooked_211.h5";
      //String fName="/work/clas12/jnp/clas_cooked_005197_irec.h5";

      // String fName="/work/clas12/jnp/irec_outfile.h5";

      // String fName="/work/clas12/jnp/instarec/irec_005197.evio.00011.h5";
      String fName="w.h5";

      String endName="_elPsup2_th0p025_matchTrack_ecalHTCC_vzElRGACuts";//_elPsup2"; //_phiCut5 eg 175-185, 10 otherwise trackChi2l350_vzl20_6SL

      double resp_threshold=0.025; //0.075
      double beamE=10.6;
        
      OnePionValidator dp = new OnePionValidator();

      //dp.process(fName,-1,endName,resp_threshold,beamE,0,0);
      //dp.process(fName,-1,endName+"_wFid",resp_threshold,beamE,1,0);
      //dp.process(fName,-1,endName+"_wFidTight",resp_threshold,beamE,2,0);


      //Fill by hand unfortunately after twig fits
      /*BarChartBuilder b = new BarChartBuilder();
      b.addEntry("L1 Trigger ",125.98,217);
      b.addEntry("L1 Trigger & Online e^- (94% / 67%)",117.882,146.336);
      b.addEntry("L1 Trigger & Offline e^- (91% / 60%)", 114.969,132.577);
      b.setTitleY("Counts");
      b.setColors(new int[]{1,2,5});
      b.setLabels(new String[]{"With Fiducial Cuts","Without Fiducial Cuts"});
      DataGroup b2 = b.build();
      TDirectory.export("plots/OnePion"+endName+".twig","/ai/validation_barchart",b2);
      TGCanvas c = new TGCanvas(1000,1000);
      //for(DataSet ds : group.getData()) c.draw(ds, "same");
      c.view().region().draw(b2);//.showLegend(0.05, 0.95);
      c.view().region().showLegend(0.05, 0.95);
      c.repaint();*/

      //Fill by hand unfortunately after twig fits
      //BarChartBuilder b = new BarChartBuilder();
      //b.addEntry("L1 Trigger & Online e^- (97% / 94% / 68%)",88.74,118.68,147.025);
      //b.addEntry("L1 Trigger & Offline e^- (96% / 91% / 61%)",86.785,114.702,132.448);
      //b.addEntry("L1 Trigger ",90.324,126.42,216.897);
      //b.setTitleY("Counts");
      //b.setColors(new int[]{2,5,1});
      //b.setLabels(new String[]{"Tight Fiducial Cuts","Loose Fiducial Cuts","No Fiducial Cuts"});
      //DataGroup b2 = b.build();
      //TDirectory.export("plots_rga/OnePion"+endName+".twig","/ai/validation_barchart",b2);
      //TGCanvas c = new TGCanvas(1000,1000);
      ////for(DataSet ds : group.getData()) c.draw(ds, "same");
      //c.view().region().draw(b2);//.showLegend(0.05, 0.95);
      //c.view().region().showLegend(0.05, 0.95);
      //c.repaint();
      
        
    }

}