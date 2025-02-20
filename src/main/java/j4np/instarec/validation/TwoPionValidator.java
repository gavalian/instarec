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

/**
 *
 * @author tyson
 */
public class TwoPionValidator {
    
    public TwoPionValidator(){
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

    public void fillRECPart(Bank RECPart, Bank trackBank, int pindex, double[] part){
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
      part[9]=getTrackSector(pindex,trackBank);
      part[10]=0;
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
          isMatched=211;
          if(pid==11 && charge==-1){ //resp>th
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

    public int getTrackSector(int pindext, Bank trackBank){
      int s=0;
      for(int i=0;i<trackBank.getRows();i++){
        int sector=trackBank.getInt("sector",i);
        int pindex=trackBank.getInt("pindex",i);
        if(pindex==pindext){
          s=sector;
        }
      }
      return s;
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


    public void process(String file, int limEvs,String endName, double threshold,double beamE, Boolean reqFids, int desired_sector){
      
      double lim_p_res=0.2; //percentage

      //output training sample size per sector
      int count = 0, nEl=0, nNotEl=0, nAllInstaEl=0;

      HipoReader r = new HipoReader();
      Event ev = new Event();
      r.open(file);
      

      /*for priting */
      Bank recpart = r.getBank("REC::Particle");
      Bank reccal = r.getBank("REC::Calorimeter");
      Bank rectrack = r.getBank("REC::Track");
      Bank triggerbank = r.getBank("RUN::config");

      int nVarsPart=11,nVarsExc=4;
      // pid, P, Theta, Phi
      double[] pim=new double[nVarsPart];
      double[] pip=new double[nVarsPart];
      double[] el=new double[nVarsPart];
      double[] exc=new double[nVarsExc];
      double[] elEs=new double[9];
      double[] elLs=new double[9];
      double[] pimEs=new double[9];
      double[] pimLs=new double[9];
      double[] pipEs=new double[9];
      double[] pipLs=new double[9];

      double targMass=getM(2212);

      String instapidelstring="";
      String instapidelfilestring="";
      if(file.contains("211")){
        instapidelstring="!";
        instapidelfilestring="_noInstaEl";
      }

      H1F hNCal = new H1F("Calo: Offline e- PID & "+instapidelstring+" Online e- PID", 6,-1.5,4.5);
      hNCal.attr().setLineColor(5);
      hNCal.attr().setLineWidth(3);
      hNCal.attr().setTitleX("Nb of Calorimeter Layers Hit");

      H1F hNCal_w = new H1F("Calo: Offline ! e- PID & "+instapidelstring+" Online e- PID", 6,-1.5,4.5);
      hNCal_w.attr().setLineColor(5);
      hNCal_w.attr().setLineWidth(3);
      hNCal_w.attr().setTitleX("Nb of Calorimeter Layers Hit");

      H1F hNCal_allInstaEl = new H1F("Calo: "+instapidelstring+" Online e- PID", 6,-1.5,4.5);
      hNCal_allInstaEl.attr().setLineColor(5);
      hNCal_allInstaEl.attr().setLineWidth(3);
      hNCal_allInstaEl.attr().setTitleX("Nb of Calorimeter Layers Hit");

      H1F hPCAL_allInstaEl = new H1F("PCAL: "+instapidelstring+" Online e- PID", 4,-1.5,2.5);
      hPCAL_allInstaEl.attr().setLineColor(5);
      hPCAL_allInstaEl.attr().setLineWidth(3);
      hPCAL_allInstaEl.attr().setTitleX("PCAL is Hit");

      H1F hECIN_allInstaEl = new H1F("ECIN: "+instapidelstring+" Online e- PID", 4,-1.5,2.5);
      hECIN_allInstaEl.attr().setLineColor(5);
      hECIN_allInstaEl.attr().setLineWidth(3);
      hECIN_allInstaEl.attr().setTitleX("ECIN is Hit");

      H1F hECOUT_allInstaEl = new H1F("ECOUT: "+instapidelstring+" Online e- PID", 4,-1.5,2.5);
      hECOUT_allInstaEl.attr().setLineColor(5);
      hECOUT_allInstaEl.attr().setLineWidth(3);
      hECOUT_allInstaEl.attr().setTitleX("ECOUT is Hit");

      H1F hHTCC_allInstaEl = new H1F("HTCC: "+instapidelstring+" Online e- PID", 4,-1.5,2.5);
      hHTCC_allInstaEl.attr().setLineColor(5);
      hHTCC_allInstaEl.attr().setLineWidth(3);
      hHTCC_allInstaEl.attr().setTitleX("HTCC is Hit");

      H1F hMM = new H1F("Mx: Offline e- PID & "+instapidelstring+" Online e- PID", 100,0,3.5);
      hMM.attr().setLineColor(5);
      hMM.attr().setLineWidth(3);
      hMM.attr().setTitleX("Mx(e'#pi^+#pi^-) [GeV]");

      H1F hMM_w = new H1F("Mx: Offline ! e- PID & "+instapidelstring+" Online e- PID", 100,0,3.5);
      hMM_w.attr().setLineColor(1);
      hMM_w.attr().setLineWidth(3);
      hMM_w.attr().setTitleX("Mx(e'#pi^+#pi^-) [GeV]");

      H1F hMM_allInstaEl = new H1F("Mx: "+instapidelstring+" Online e- PID", 100,0,3.5);
      hMM_allInstaEl.attr().setLineColor(2);
      hMM_allInstaEl.attr().setLineWidth(3);
      hMM_allInstaEl.attr().setTitleX("Mx(e'#pi^+#pi^-) [GeV]");

      H1F hIM = new H1F("M: Offline e- PID & "+instapidelstring+" Online e- PID", 100,0,3.5);
      hIM.attr().setLineColor(5);
      hIM.attr().setLineWidth(3);
      hIM.attr().setTitleX("M(#pi^+#pi^-) [GeV]");

      H1F hIM_w = new H1F("M: Offline ! e- PID & "+instapidelstring+" Online e- PID", 100,0,3.5);
      hIM_w.attr().setLineColor(1);
      hIM_w.attr().setLineWidth(3);
      hIM_w.attr().setTitleX("M(#pi^+#pi^-) [GeV]");

      H1F hIM_allInstaEl = new H1F("M: "+instapidelstring+" Online e- PID", 100,0,3.5);
      hIM_allInstaEl.attr().setLineColor(2);
      hIM_allInstaEl.attr().setLineWidth(3);
      hIM_allInstaEl.attr().setTitleX("M(#pi^+#pi^-) [GeV]");

      H1F hPDif = new H1F("e^- P - #pi^- P "+instapidelstring+" Online e- PID", 100,-5,5);
      hPDif.attr().setLineColor(2);
      hPDif.attr().setLineWidth(3);
      hPDif.attr().setTitleX("e^- P - #pi^- P [GeV]");

      H1F hPDif_all = new H1F("e^- P - #pi^- P "+instapidelstring+" Online e- PID (all phase space)", 100,-5,5);
      hPDif_all.attr().setLineColor(2);
      hPDif_all.attr().setLineWidth(3);
      hPDif_all.attr().setTitleX("e^- P - #pi^- P [GeV]");

      H2F hpimPelP = new H2F("e^- P vs #pi^- P "+instapidelstring+" Online e- PID",50,0,10,50,0,10);
      hpimPelP.attr().setTitleX("#pi^- P [GeV]");
      hpimPelP.attr().setTitleY("e^- P [GeV]");

      H2F helPTheta = new H2F("e^- #theta P: Offline e- PID & "+instapidelstring+" Online e- PID",50,0,10,50,0,40);
      helPTheta.attr().setTitleX("P [GeV]");
      helPTheta.attr().setTitleY("#theta [Deg]");

      H2F helPTheta_w = new H2F("e^- #theta P: Offline ! e- PID & "+instapidelstring+" Online e- PID",50,0,10,50,0,40);
      helPTheta_w.attr().setTitleX("P [GeV]");
      helPTheta_w.attr().setTitleY("#theta [Deg]");

      H2F hpimPTheta_allInstaEl = new H2F("#pi^- #theta P:  "+instapidelstring+" Online e- PID ",50,0,10,50,0,40);
      hpimPTheta_allInstaEl.attr().setTitleX("P [GeV]");
      hpimPTheta_allInstaEl.attr().setTitleY("#theta [Deg]");

      H2F hpipPTheta_allInstaEl = new H2F("#pi^+ #theta P:  "+instapidelstring+" Online e- PID ",50,0,10,50,0,40);
      hpipPTheta_allInstaEl.attr().setTitleX("P [GeV]");
      hpipPTheta_allInstaEl.attr().setTitleY("#theta [Deg]");


      H2F helPTheta_allInstaEl = new H2F("e^- #theta P:  "+instapidelstring+" Online e- PID",50,0,10,50,0,40);
      helPTheta_allInstaEl.attr().setTitleX("P [GeV]");
      helPTheta_allInstaEl.attr().setTitleY("#theta [Deg]");

      H2F helPTheta_allInstaEl_all = new H2F("e^- #theta P:  "+instapidelstring+" Online e- PID (all phase space)",50,0,10,50,0,40);
      helPTheta_allInstaEl_all.attr().setTitleX("P [GeV]");
      helPTheta_allInstaEl_all.attr().setTitleY("#theta [Deg]");

      H2F hpimPTheta_allInstaEl_all = new H2F("#pi^- #theta P:  "+instapidelstring+" Online e- PID (all phase space)",50,0,10,50,0,40);
      hpimPTheta_allInstaEl_all.attr().setTitleX("P [GeV]");
      hpimPTheta_allInstaEl_all.attr().setTitleY("#theta [Deg]");

      H2F hpipPTheta_allInstaEl_all = new H2F("#pi^+ #theta P:  "+instapidelstring+" Online e- PID (all phase space)",50,0,10,50,0,40);
      hpipPTheta_allInstaEl_all.attr().setTitleX("P [GeV]");
      hpipPTheta_allInstaEl_all.attr().setTitleY("#theta [Deg]");

      H2F helXY = new H2F("e^- Y vs X: Offline e- PID & "+instapidelstring+" Online e- PID",100,-250,250,100,-250,250);
      helXY.attr().setTitleX("X [cm]");
      helXY.attr().setTitleY("Y [cm]");

      H2F helXY_w = new H2F("e^- Y vs X: Offline ! e- PID & "+instapidelstring+" Online e- PID",100,-250,250,100,-250,250);
      helXY_w.attr().setTitleX("X [cm]");
      helXY_w.attr().setTitleY("Y [cm]");

      H2F helXY_allInstaEl = new H2F("e^- Y vs X:  "+instapidelstring+" Online e- PID",100,-250,250,100,-250,250);
      helXY_allInstaEl.attr().setTitleX("X [cm]");
      helXY_allInstaEl.attr().setTitleY("Y [cm]");

      H2F helPSF = new H2F("e^- SF vs P: Offline e- PID & "+instapidelstring+" Online e- PID",50,0,10,50,0,0.5);
      helPSF.attr().setTitleX("P [GeV]");
      helPSF.attr().setTitleY("Sampling Fraction");

      H2F helPSF_w = new H2F("e^- SF vs P: Offline ! e- PID & "+instapidelstring+" Online e- PID",50,0,10,50,0,0.5);
      helPSF_w.attr().setTitleX("P [GeV]");
      helPSF_w.attr().setTitleY("Sampling Fraction");

      H2F helPSF_allInstaEl = new H2F("e^- SF vs P:  "+instapidelstring+" Online e- PID",50,0,10,50,0,0.5);
      helPSF_allInstaEl.attr().setTitleX("P [GeV]");
      helPSF_allInstaEl.attr().setTitleY("Sampling Fraction");

      H2F helLVSF = new H2F("e^- SF vs LV: Offline e- PID & "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
      helLVSF.attr().setTitleX("LV [cm]");
      helLVSF.attr().setTitleY("Sampling Fraction");

      H2F helLVSF_w = new H2F("e^- SF vs LV: Offline ! e- PID & "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
      helLVSF_w.attr().setTitleX("LV [cm]");
      helLVSF_w.attr().setTitleY("Sampling Fraction");

      H2F helLVSF_allInstaEl = new H2F("e^- SF vs LV:  "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
      helLVSF_allInstaEl.attr().setTitleX("LV [cm]");
      helLVSF_allInstaEl.attr().setTitleY("Sampling Fraction");

      H2F helLWSF = new H2F("e^- SF vs LW: Offline e- PID & "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
      helLWSF.attr().setTitleX("LW [cm]");
      helLWSF.attr().setTitleY("Sampling Fraction");

      H2F helLWSF_w = new H2F("e^- SF vs LW: Offline ! e- PID & "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
      helLWSF_w.attr().setTitleX("LW [cm]");
      helLWSF_w.attr().setTitleY("Sampling Fraction");

      H2F helLWSF_allInstaEl = new H2F("e^- SF vs LW:  "+instapidelstring+" Online e- PID",50,0,100,50,0,0.5);
      helLWSF_allInstaEl.attr().setTitleX("LW [cm]");
      helLWSF_allInstaEl.attr().setTitleY("Sampling Fraction");

      //don't count if set limEvs=-1
      if(limEvs==-1){count=-2;}
      
      while(r.hasNext() && count<limEvs){
        //don't count if set limEvs=-1
        if(limEvs!=-1){
          count++;
        }

        r.nextEvent(ev);
        Leaf pred_part = ev.readLeaf(1,12,32200,1);
        ev.read(recpart);
        ev.read(reccal);
        ev.read(rectrack);
        ev.read(triggerbank);

        int hasTrig=hasTriggerEl(triggerbank);
        

        //System.out.println("\n\nNew Event");
        //pred_part.print();

        for(int row=0;row<recpart.getRows();row++){
          for(int row2=(row+1);row2<recpart.getRows();row2++){
            for(int row3=0;row3<recpart.getRows();row3++){

              cleanArr(pim,nVarsPart);
              cleanArr(pip,nVarsPart);
              cleanArr(el,nVarsPart);
              cleanArr(exc,nVarsExc);
              cleanArr(elEs,9);
              cleanArr(elLs,9);
              cleanArr(pimEs,9);
              cleanArr(pimLs,9);
              cleanArr(pipEs,9);
              cleanArr(pipLs,9);
              fillRECPart(recpart,rectrack,row,el);
              getCalInfo(reccal, row, el[1], elEs, elLs);
              int matchel=matchTracks(pred_part,el,lim_p_res,threshold);
              fillRECPart(recpart,rectrack,row2,pim);
              getCalInfo(reccal, row2, pim[1], pimEs, pimLs);
              int matchpim=matchTracks(pred_part,pim,lim_p_res,threshold);
              fillRECPart(recpart,rectrack,row3,pip);
              getCalInfo(reccal, row3, pip[1], pipEs, pipLs);
              int matchpip=matchTracks(pred_part,pip,lim_p_res,threshold);

              int desired_sector_e=desired_sector;
              int sector=(int)el[9];
              //System.out.printf("p pid %f px %f py %f pz %f status %f sector %f charge %f match %d\n",part[0],part[6],part[7],part[8],part[4],part[9],part[5],match);

              //use all sectors if desired sector is 0
              if(desired_sector==0){
                desired_sector_e=sector;
              }
            
              Boolean fid=passFid(elLs);
              if(!reqFids){
                fid=true;
              }
            
              int hasEl=0,hasPim=0,hasPip=0;
              if(file.contains("211")){
                //only one pos pion
                if(pip[5]==1 && matchpip!=-1 && pip[0]==211){
                  hasPip=1;
                }
                if(el[5]==-1 && matchel!=-1 && desired_sector_e==sector && fid){
                  if(el[0]==-211 || el[0]==11){
                    hasEl=1;
                  }
                }
                if(pim[5]==-1 && matchpim!=-1 && pim[0]==-211){
                  hasPim=1;
                }
              } else{
                if(pip[5]==1 && matchpip!=-1 && pip[0]==211){
                  hasPip=1;
                }
                if(el[5]==-1 && matchel==11 && desired_sector_e==sector && fid){
                  if(el[0]==-211 || el[0]==11){
                    hasEl=1;
                  }
                }
                if(pim[5]==-1 && matchpim!=-1 && matchpim!=11 && pim[0]==-211){
                    hasPim=1;
                }
              }//file 211 or not

              if(hasEl==1 && hasPim==1 && hasPip==1 && hasTrig==1){
                calcExc(pip,pim, el, exc,beamE);
      
                if(exc[1]>0.2 && el[1]>2){
                  //System.out.println("\n\nFound particles");
                  //System.out.printf("el pid %f px %f py %f pz %f status %f sector %f charge %f\n",el[0],el[6],el[7],el[8],el[4],el[9],el[5]);
                  //System.out.printf("pi- pid %f px %f py %f pz %f status %f sector %f charge %f\n",pim[0],pim[6],pim[7],pim[8],pim[4],pim[9],pim[5]);
                  //System.out.printf("pi+ pid %f px %f py %f pz %f status %f sector %f charge %f\n",pip[0],pip[6],pip[7],pip[8],pip[4],pip[9],pip[5]);
      
                  
                  int elNCal=0, hasPCAL=0,hasECIN=0,hasECOUT=0,hasHTCC=0;
                  if(elEs[0]>0.01){elNCal++;hasPCAL=1;}
                  if(elEs[1]>0.01){elNCal++;hasECIN=1;}
                  if(elEs[2]>0.01){elNCal++;hasECOUT=1;}
                  if(el[10]>0.0){hasHTCC=1;}
                  
      
                  nAllInstaEl++;
                  hIM_allInstaEl.fill(exc[0]);
                  hMM_allInstaEl.fill(exc[1]);
                  hpimPelP.fill(pim[1],el[1]);
                  hPDif_all.fill(el[1]-pim[1]);
                  //if(exc[0]>0.7 && exc[0]<0.86){
                  helPTheta_allInstaEl_all.fill(el[1],el[2]);
                  hpimPTheta_allInstaEl_all.fill(pim[1],pim[2]);
                  hpipPTheta_allInstaEl_all.fill(pip[1],pip[2]);
                  if(exc[1]>0.85 && exc[1]<1){
                    hPDif.fill(el[1]-pim[1]);
                    helPTheta_allInstaEl.fill(el[1],el[2]);
                    helPSF_allInstaEl.fill(el[1],elEs[3]);
                    helXY_allInstaEl.fill(elEs[7],elEs[8]);
                    hpimPTheta_allInstaEl.fill(pim[1],pim[2]);
                    hpipPTheta_allInstaEl.fill(pip[1],pip[2]);
                    helLVSF_allInstaEl.fill(elLs[1],elEs[3]);
                    helLWSF_allInstaEl.fill(elLs[2],elEs[3]);
                    hNCal_allInstaEl.fill(elNCal);
                    hPCAL_allInstaEl.fill(hasPCAL);
                    hECIN_allInstaEl.fill(hasECIN);
                    hECOUT_allInstaEl.fill(hasECOUT);
                    hHTCC_allInstaEl.fill(hasHTCC);
                  }
                  if(el[0]==11){
                    hIM.fill(exc[0]);
                    hMM.fill(exc[1]);
                    //if(exc[0]>0.7 && exc[0]<0.86){
                    if(exc[1]>0.85 && exc[1]<1){
                      helPTheta.fill(el[1],el[2]);
                      helPSF.fill(el[1],elEs[3]);
                      helXY.fill(elEs[7],elEs[8]);
                      helLVSF.fill(elLs[1],elEs[3]);
                      helLWSF.fill(elLs[2],elEs[3]);
                      hNCal.fill(elNCal);
                    }
                      nEl++;
                  } else {
                    hIM_w.fill(exc[0]);
                    hMM_w.fill(exc[1]);
                    //if(exc[0]>0.7 && exc[0]<0.86){
                    if(exc[1]>0.85 && exc[1]<1){
                      helPTheta_w.fill(el[1],el[2]);
                      helPSF_w.fill(el[1],elEs[3]);
                      helXY_w.fill(elEs[7],elEs[8]);
                      helLVSF_w.fill(elLs[1],elEs[3]);
                      helLWSF_w.fill(elLs[2],elEs[3]);
                      hNCal_w.fill(elNCal);
                    }
                      nNotEl++;
                  }
                }
              }

            }//loop over rec parts 3
          }// loop over rec parts 2
        }//loop over rec parts 1

      }//read file

      if(desired_sector==0){
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hMM);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hMM_w);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hMM_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hIM);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hIM_w);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hIM_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hPDif);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hPDif_all);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hpimPelP);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helPTheta);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helPSF);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helLVSF);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helLWSF);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helXY);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hNCal);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helPTheta_w);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helPSF_w);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helLVSF_w);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helLWSF_w);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helXY_w);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hNCal_w);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helPTheta_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hpimPTheta_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hpipPTheta_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helPSF_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helLVSF_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helLWSF_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helXY_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,helPTheta_allInstaEl_all);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hpimPTheta_allInstaEl_all);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hpipPTheta_allInstaEl_all);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hNCal_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hPCAL_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hECIN_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hECOUT_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation"+instapidelfilestring,hHTCC_allInstaEl);
      } else {
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hMM);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hMM_w);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hMM_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hIM);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hIM_w);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hIM_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hPDif);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hPDif_all);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hpimPelP);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helPTheta);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helPSF);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helLVSF);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helLWSF);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helXY);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hNCal);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helPTheta_w);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helPSF_w);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helLVSF_w);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helLWSF_w);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helXY_w);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hNCal_w);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helPTheta_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hpimPTheta_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hpipPTheta_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helPSF_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helLVSF_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helLWSF_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helXY_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hNCal_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,helPTheta_allInstaEl_all);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hpimPTheta_allInstaEl_all);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hpipPTheta_allInstaEl_all);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hPCAL_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hECIN_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hECOUT_allInstaEl);
        TDirectory.export("plots/TwoPion"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,hHTCC_allInstaEl);
      }

      //System.out.printf("Number of Events with Insta And Rec El %d, Not Rec el %d \n\n",nEl,nNotEl);
      
    }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar
    //read plots in j4shell with eg TwigStudio.browser("plots/TwoPion.twig")
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting elastic scattering validator ");
      String fName_11="/Users/tyson/data_repo/trigger_data/rga/instarec_skims/cooked_11.h5";
      String fName_211="/Users/tyson/data_repo/trigger_data/rga/instarec_skims/cooked_211.h5";

      String endName="_elPsup2"; //_phiCut5 eg 175-185, 10 otherwise

      double resp_threshold=0.075;
      double beamE=10.6;
        
      TwoPionValidator dp = new TwoPionValidator();
      dp.process(fName_11,-1,endName,resp_threshold,beamE,false,0);
      dp.process(fName_11,-1,endName+"_wFid",resp_threshold,beamE,true,0);
      dp.process(fName_211,-1,endName,resp_threshold,beamE,false,0);
      dp.process(fName_211,-1,endName+"_wFid",resp_threshold,beamE,true,0);


      //Fill by hand unfortunately after twig fits
      /*BarChartBuilder b = new BarChartBuilder();
      b.addEntry("L1 e- ",93.,154.);
      b.addEntry("L1 e- & Online e- (100% / 78%)",93,120);
      b.addEntry("L1 e- & Online e- & Offline e- (92% / 65%)", 86,100);
      b.setTitleY("Counts");
      b.setColors(new int[]{1,2,5});
      b.setLabels(new String[]{"With Fiducial Cuts","Without Fiducial Cuts"});
      DataGroup b2 = b.build();
      TDirectory.export("plots/TwoPion"+endName+".twig","/ai/validation_barchart",b2);
      TGCanvas c = new TGCanvas(1000,1000);
      //for(DataSet ds : group.getData()) c.draw(ds, "same");
      c.view().region().draw(b2);//.showLegend(0.05, 0.95);
      c.view().region().showLegend(0.05, 0.95);
      c.repaint();*/
      
        
    }

}