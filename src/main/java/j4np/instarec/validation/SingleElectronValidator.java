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
public class SingleElectronValidator {
    
    public SingleElectronValidator(){
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


    public void process(String file, int limEvs,String endName, double threshold, int reqFids, int desired_sector){
      
      double lim_p_res=0.2; //percentage

      //output training sample size per sector
      int count = 0, nEl=0, nNotEl=0, nAllInstaEl=0;

      HipoReader r = new HipoReader();
      Event ev = new Event();
      //r.setTags(2211); //2211 for two pion tagger
      r.open(file);
      

      /*for printing */
      Bank recpart = r.getBank("REC::Particle");
      Bank reccal = r.getBank("REC::Calorimeter");
      Bank rectrack = r.getBank("REC::Track");
      Bank triggerbank = r.getBank("RUN::config");

      Leaf pred_part = new Leaf(32, 3, "i", 1200);

      int nVarsPart=14;
      // pid, P, Theta, Phi=
      double[] el=new double[nVarsPart];
      double[] elEs=new double[9];
      double[] elLs=new double[9];

      double targMass=getM(2212);

      Vector<H1F> histos = new Vector<>();
      Map<String, Integer> names = new HashMap<>();

      Vector<H2F> histos2D = new Vector<>();
      Map<String, Integer> names2D = new HashMap<>();

      int nHistos=0,nHistos2D=0;
      for(int i=0;i<3;i++){
        String instapidelstring="All";
        String instapidelfilestring="_All";
        int col=1;
        if(i==1){
          instapidelstring="wElID";
          instapidelfilestring="_wElID";
          col=2;
        } else if(i==2){
          instapidelstring="notElID";
          instapidelfilestring="_notElID";
          col=5;
        }

        H2F helPTheta = new H2F("e^- #theta P: "+instapidelstring,50,0,10,50,0,40);
        helPTheta.attr().setTitleX("P [GeV]");
        helPTheta.attr().setTitleY("#theta [Deg]");
        names2D.put("helPTheta"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helPTheta);

        H2F helSFLV = new H2F("e^- Sampling Fraction vs LV "+instapidelstring,50,0,50,50,0,0.5);
        helSFLV.attr().setTitleX("LV [cm]");
        helSFLV.attr().setTitleY("Sampling Fraction");
        names2D.put("helSFLV"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helSFLV);

        H2F helSFLW = new H2F("e^- Sampling Fraction vs LW "+instapidelstring,50,0,50,50,0,0.5);
        helSFLW.attr().setTitleX("LW [cm]");
        helSFLW.attr().setTitleY("Sampling Fraction");
        names2D.put("helSFLW"+instapidelstring,nHistos2D);
        nHistos2D++;
        histos2D.add(helSFLW);

        H1F helP = new H1F("e^- P "+instapidelstring, 100,0,10);
        helP.attr().setLineColor(col);
        helP.attr().setLineWidth(3);
        helP.attr().setTitleX("e^- P [GeV]");
        names.put("helP"+instapidelstring,nHistos);
        nHistos++;
        histos.add(helP);
      }

      

      //don't count if set limEvs=-1
      if(limEvs==-1){count=-2;}
      
      while(r.hasNext() && count<limEvs){
        

        r.nextEvent(ev);
        //Leaf pred_part = ev.readLeaf(1,12,32,3);
        ev.read(pred_part,32,3);
        ev.read(recpart);
        ev.read(reccal);
        ev.read(rectrack);
        ev.read(triggerbank);

        int hasTrig=hasTriggerEl(triggerbank);
        

        //System.out.printf("\n\nNew Event, %d \n",pred_part.getSize());
        //pred_part.print();
        //recpart.show();
        //reccal.show();

        if(pred_part.getSize()>0){

          for(int row=0;row<recpart.getRows();row++){

            cleanArr(el,nVarsPart);
            cleanArr(elEs,9);
            cleanArr(elLs,9);

            fillRECPart(recpart,rectrack,row,el);
            int matchel=0;

            //System.out.printf("p pid %f px %f py %f pz %f status %f sector %f charge %f match %d\n",part[0],part[6],part[7],part[8],part[4],part[9],part[5],match);

            int desired_sector_el=desired_sector;
            int sector=(int)el[9];
            //use all sectors if desired sector is 0
            if(desired_sector==0){
              desired_sector_el=sector;
            }
                
            int hasElCandi=0;
            if(el[5]==-1 && el[0]==11 ){ //&& matchel!=-1 el[11] is track chi^2 && el[11]<350 && Math.abs(el[13])<20 && el[12]==6
              matchel=matchTracks(pred_part,el,lim_p_res,threshold);
              if(matchel!=-1){
                hasElCandi=1;
                getCalInfo(reccal, row, el[1], elEs, elLs);
              }
              
            }

            /*if(hasElCandi==1){
              System.out.println("\n\nFound particle");
              System.out.printf("el pid %f px %f py %f pz %f status %f sector %f charge %f\n",el[0],el[6],el[7],el[8],el[4],el[9],el[5]);
            }*/

            if(hasElCandi==1 && hasTrig==1){
              if(limEvs!=-1){
                count++;
              }
              //System.out.println("\n\nFound particles");
              //System.out.printf("el pid %f px %f py %f pz %f status %f sector %f charge %f\n",el[0],el[6],el[7],el[8],el[4],el[9],el[5]);
              //System.out.printf("pi- pid %f px %f py %f pz %f status %f sector %f charge %f\n",pim[0],pim[6],pim[7],pim[8],pim[4],pim[9],pim[5]);
              //System.out.printf("pi+ pid %f px %f py %f pz %f status %f sector %f charge %f\n",pip[0],pip[6],pip[7],pip[8],pip[4],pip[9],pip[5]);

              for(int j=0;j<3;j++){
                String instapidelstring="All";
                int cont=1;
                if(j==1){
                  instapidelstring="wElID";
                  if(matchel==11){cont=1;}
                  else{cont=0;}
                } else if (j==2){
                  instapidelstring="notElID";
                  if(matchel!=11){cont=1;}
                  else{cont=0;}
                }
                if(cont==1){
                  fillHisto2D(names2D,histos2D,"helPTheta"+instapidelstring,el[1],el[2]);
                  fillHisto2D(names2D,histos2D,"helSFLV"+instapidelstring,elLs[1],elEs[3]);
                  fillHisto2D(names2D,histos2D,"helSFLW"+instapidelstring,elLs[2],elEs[3]);
                  fillHisto(names,histos,"helP"+instapidelstring,el[1]);          
                }      
              }//fill histos 
            }//has candi
          }//loop over rec parts 1
        }//check bank not empty

      }//read file
  
      for (Map.Entry<String, Integer> entry : names.entrySet()) {
        String name=entry.getKey();
        int index=entry.getValue();

        String instapidelfilestring="";
        if(name.contains("All")){instapidelfilestring="_All";}
        if(name.contains("wElID")){instapidelfilestring="_wElID";}
        if(name.contains("notElID")){instapidelfilestring="_notElID";}

        if(desired_sector==0){
          TDirectory.export("plots_rga/SingleElectron"+endName+".twig","/ai/validation"+instapidelfilestring,histos.get(index));
        } else {
          TDirectory.export("plots_rga/SingleElectron"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,histos.get(index));
        }
      }

      for (Map.Entry<String, Integer> entry : names2D.entrySet()) {
        String name=entry.getKey();
        int index=entry.getValue();

        String instapidelfilestring="";
        if(name.contains("All")){instapidelfilestring="_All";}
        if(name.contains("wElID")){instapidelfilestring="_wElID";}
        if(name.contains("notElID")){instapidelfilestring="_notElID";}

        if(desired_sector==0){
          TDirectory.export("plots_rga/SingleElectron"+endName+".twig","/ai/validation"+instapidelfilestring,histos2D.get(index));
        } else {
          TDirectory.export("plots_rga/SingleElectron"+endName+"_Sector"+String.valueOf(desired_sector)+".twig","/ai/validation"+instapidelfilestring,histos2D.get(index));
        }
      }
      
    }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar
    //read plots in j4shell with eg TwigStudio.browser("plots/SingleElectron.twig")
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting Single Electron validator ");
      //String fName_11="/Users/tyson/data_repo/trigger_data/rga/instarec_skims/cooked_11.h5";
      //String fName_211="/Users/tyson/data_repo/trigger_data/rga/instarec_skims/cooked_211.h5";
      //String fName="/work/clas12/jnp/clas_cooked_005197_irec.h5";

      // String fName="/work/clas12/jnp/irec_outfile.h5";

      // String fName="/work/clas12/jnp/instarec/irec_005197.evio.00011.h5";
      String fName="wvalid.h5";

      String endName="_noPimID_matchTrack";//_elPsup2"; //_phiCut5 eg 175-185, 10 otherwise trackChi2l350_vzl20_6SL

      double resp_threshold=0.025; //0.075
      double beamE=10.6;
        
      SingleElectronValidator dp = new SingleElectronValidator();

      dp.process(fName,200000,endName,resp_threshold,0,0);
      //dp.process(fName,200000,endName+"_wFid",resp_threshold,1,0);
      //dp.process(fName,200000,endName+"_wFidTight",resp_threshold,2,0);


      fName="wRadPhotonsSmall.h5";
      endName="_radPhotonTraining"+endName;

      //dp.process(fName,200000,endName,resp_threshold,0,0);
      //dp.process(fName,200000,endName+"_wFid",resp_threshold,1,0);
      //dp.process(fName,200000,endName+"_wFidTight",resp_threshold,2,0);
      
        
    }

}