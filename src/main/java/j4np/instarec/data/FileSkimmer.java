/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.data;

import j4np.hipo5.data.Bank;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.hipo5.io.HipoWriter;
import java.io.File;


/**
 *
 * @author tyson
 */
public class FileSkimmer {

    
    
    public FileSkimmer(){
    }

    public Boolean radiatedPhoton(Bank RECPart, int pindex){

      double pz = RECPart.getFloat("pz", pindex);
      double px = RECPart.getFloat("px", pindex);
      double py = RECPart.getFloat("py", pindex);
      double p=Math.sqrt(px*px+py*py+pz*pz);
      double Theta = Math.acos(pz / p)*(180/Math.PI);// Math.atan2(Math.sqrt(px*px+py*py),pz);
      double Phi = Math.atan2(py, px)*(180/Math.PI);
      for (short i = 0; i < RECPart.getRows(); i++) {
        int charge = RECPart.getByte("charge", i);
        double pzphoton = RECPart.getFloat("pz", i);
        double pxphoton = RECPart.getFloat("px", i);
        double pyphoton = RECPart.getFloat("py", i);
        double pphoton=Math.sqrt(pxphoton*pxphoton+pyphoton*pyphoton+pzphoton*pzphoton);
        double Thetaphoton = Math.acos(pzphoton / pphoton)*(180/Math.PI);// Math.atan2(Math.sqrt(pxphoton*pxphoton+pyphoton*pyphoton),pzphoton);
        double Phiphoton = Math.atan2(pyphoton, pxphoton)*(180/Math.PI);
        //just use charge and not pid because some photons missided as neutrons
        if(Math.abs(Theta-Thetaphoton)<3 && charge==0){
          return true;
        }
      }

      return false;

    }

    public Boolean hasPhotonElectronPair(Bank RECPart){

      for(short pindex=0; pindex<RECPart.getRows();pindex++){
        int charge = RECPart.getInt("charge", pindex);
        short status = RECPart.getShort("status", pindex);
        //use charge not pid cos want to allow cases for missIDed e- by rec
        if(charge==-1 && Math.abs(status)>=2000 & Math.abs(status)<4000){
          double pz = RECPart.getFloat("pz", pindex);
          double px = RECPart.getFloat("px", pindex);
          double py = RECPart.getFloat("py", pindex);
          double p=Math.sqrt(px*px+py*py+pz*pz);
          double Theta = Math.acos(pz / p)*(180/Math.PI);// Math.atan2(Math.sqrt(px*px+py*py),pz);
          double Phi = Math.atan2(py, px)*(180/Math.PI);
          for (short i = 0; i < RECPart.getRows(); i++) {
            int chargephoton = RECPart.getInt("charge", i);
            double pzphoton = RECPart.getFloat("pz", i);
            double pxphoton = RECPart.getFloat("px", i);
            double pyphoton = RECPart.getFloat("py", i);
            double pphoton=Math.sqrt(pxphoton*pxphoton+pyphoton*pyphoton+pzphoton*pzphoton);
            double Thetaphoton = Math.acos(pzphoton / pphoton)*(180/Math.PI);// Math.atan2(Math.sqrt(pxphoton*pxphoton+pyphoton*pyphoton),pzphoton);
            double Phiphoton = Math.atan2(pyphoton, pxphoton)*(180/Math.PI);
            //just use chargephoton and not pid because some photons missided as neutrons
            if(Math.abs(Theta-Thetaphoton)<3 && chargephoton==0 && i!=pindex){
              return true;
            }
          }
        }
      }

      return false;

    }

    public int hasFDPID(Bank RECPart, int pidw){
      for (short i = 0; i < RECPart.getRows(); i++) {
        int pid = RECPart.getInt("pid", i);
        short status = RECPart.getShort("status", i);
        if(pid==pidw){
          if(Math.abs(status)>=2000 && Math.abs(status)<4000){
            return i;
          }
        }
      }
      return -1;

    }

    public int countChargeFD(Bank RECPart, int chargew){
      int countCh=0;
      for (short i = 0; i < RECPart.getRows(); i++) {
        int charge = RECPart.getInt("charge", i);
        short status = RECPart.getShort("status", i);
        if(charge==chargew){
          if(Math.abs(status)>=2000 && Math.abs(status)<4000){
            countCh++;
          }
        }
      }
      return countCh;

    }

    public void cleanArr(double[] arr, double length){
      for(int i=0;i<length;i++){arr[i]=0;}
    }

    public void copyArr(double[] arr,double[] arr_cp, double length){
      for(int i=0;i<length;i++){arr_cp[i]=arr[i];}
    }

    public double square(double a){
      return a*a;
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

    public void calcExc(double[] pair1, double[] pair2, double[] el, double[] exc,double beamE, int squareMM){
      
      double elE=Math.sqrt(square(el[1])+square(getM((int)el[0])));
      double pair1E=Math.sqrt(square(pair1[1])+square(getM((int)pair1[0])));
      double pair2E=Math.sqrt(square(pair2[1])+square(getM((int)pair2[0])));
      double pM=getM(2212);

      double IM = Math.sqrt(square(pair1E+pair2E)- ( square(pair1[6]+pair2[6]) + square(pair1[7]+pair2[7]) + square(pair1[8]+pair2[8]) ));
      double px_m = -1.0*(el[6]+pair1[6]+pair2[6]);
      double py_m = -1.0*(el[7]+pair1[7]+pair2[7]);
      double pz_m = beamE - (el[8]+pair1[8]+pair2[8]);
      double p_m=Math.sqrt(square(px_m) + square(py_m)+square(pz_m));
      double pxp_m = px_m/p_m;
      double pyp_m = py_m/p_m;
      double E_m = beamE + getM(2212) - (elE+pair1E+pair2E);
      double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));

      exc[0]=IM;
      exc[1]=MM2;
      if(squareMM==1){
        exc[1]=Math.sqrt(MM2);
      }
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
      part[9]=vz;
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

    public Boolean passFid(double[] Ls, double part[],int str){

      boolean vz=false,looseCal=false,tightCal=false;

      //for str 0:
      if(part[9]<20 && part[9]>(-13)){
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

    public void skimTwoPion(String inDir,String outFile,double beamE, double targetMass, double MMUpperLim, int cutRho, int limEvs, int FidStr){

      File folder = new File(inDir);

      File[] files = folder.listFiles();

      HipoReader rtmp = new HipoReader(inDir+files[0].getName());
      HipoWriter w = HipoWriter.create(outFile, rtmp);

      long count=0, countWrite=0, nFile=0;

      //don't count if set limEvs=-1
      if(limEvs==-1){countWrite=-2;}

      for (File file : files) {
        nFile++;
        if (file.isFile() && countWrite<limEvs) {
          HipoReader r = new HipoReader(inDir+file.getName());
          Event ev = new Event();
          Event ev_out = new Event();

          Bank recpart = r.getBank("REC::Particle");
          Bank reccal = r.getBank("REC::Calorimeter");
          Bank dc = r.getBank("DC::tdc");
          Bank ecal = r.getBank("ECAL::adc");
          Bank ftoftdc = r.getBank("FTOF::tdc");
          Bank ftofadc = r.getBank("FTOF::adc");
          Bank htcc = r.getBank("HTCC::adc");
          Bank rechtcc = r.getBank("REC::Cherenkov");
          Bank recftof = r.getBank("REC::Scintillator");
          Bank rectrack = r.getBank("REC::Track");
          Bank runconfig = r.getBank("RUN::config");

          int nVarsPart=10, nVarsExc=4;
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
          while(r.hasNext() && countWrite<limEvs){

            r.nextEvent(ev);
            ev.read(recpart);
            ev.read(reccal);
            ev.read(rechtcc);
            ev.read(recftof);
            ev.read(rectrack);
            ev.read(runconfig);
            ev.read(ecal);
            ev.read(htcc);
            ev.read(ftoftdc);
            ev.read(ftofadc);
            ev.read(dc);

            int hasTrig=hasTriggerEl(runconfig);

            int addEv=0;

            for(int row=0;row<recpart.getRows();row++){
              for(int row2=(row+1);row2<recpart.getRows();row2++){
                for(int row3=0;row3<recpart.getRows();row3++){

                  cleanArr(pim,nVarsPart);
                  cleanArr(pip,nVarsPart);
                  cleanArr(el,nVarsPart);
                  cleanArr(exc,nVarsExc);
                  

                  fillRECPart(recpart,rectrack,row,el);
                  fillRECPart(recpart,rectrack,row2,pim);
                  fillRECPart(recpart,rectrack,row3,pip);
                  int hasElCandi=0,hasPim=0,hasPip=0;

                  //use el fid to get good sample of negative pions
                  Boolean fid=true;
                  
                  //fids adds extra detector, makes things slow
                  //ignore with fidstrength=-1
                  if(FidStr!=-1){
                    cleanArr(elEs,9);
                    cleanArr(elLs,9);
                    cleanArr(pimEs,9);
                    cleanArr(pimLs,9);
                    cleanArr(pipEs,9);
                    cleanArr(pipLs,9);
                    getCalInfo(reccal, row, el[1], elEs, elLs);
                    getCalInfo(reccal, row2, pim[1], pimEs, pimLs);
                    getCalInfo(reccal, row3, pip[1], pipEs, pipLs);
                    passFid(elLs,el,FidStr);
                  }
                  

                  if(pip[5]==1 && pip[0]==211){
                    hasPip=1;
                  }
                  if(pim[5]==-1 ){ //&& pim[0]==-211
                    hasPim=1;
                  }

                  if(el[5]==-1 && fid){ //&& matchel!=-1 el[11] is track chi^2 && el[11]<350 && Math.abs(el[13])<20 && el[12]==6
                    hasElCandi=1;
                  }

                  if(hasElCandi==1 && hasPim==1 && hasPip==1 && hasTrig==1){
                    calcExc(pip,pim, el, exc,beamE,1);
                    //System.out.printf("\nMM %f",exc[1]);
                    if(cutRho==0){
                      if(exc[1]>0.2 && exc[1]<MMUpperLim){
                        //System.out.printf(" after MM %f \n",exc[1]);
                        addEv=1;
                      }
                    } else{
                      if(exc[1]>0.2 && exc[1]<MMUpperLim){
                        if(exc[0]>0.6 && exc[0]<0.9){
                          //System.out.printf(" after MM %f \n",exc[1]);
                          addEv=1;
                        }
                      }
                    }//cut on rho or not
                    
                  }//have good candidates

                }//read pi+
              }//read pi-
            }//read el

            if(addEv==1){
              ev_out.reset();
              ev_out.write(recpart);
              ev_out.write(reccal);
              ev_out.write(rechtcc);
              ev_out.write(recftof);
              ev_out.write(rectrack);
              ev_out.write(runconfig);
              ev_out.write(ecal);
              ev_out.write(htcc);
              ev_out.write(ftoftdc);
              ev_out.write(ftofadc);
              ev_out.write(dc);
              w.addEvent(ev_out);
              count++;
              if(limEvs!=-1){
                countWrite++;
              }
            }

          }//loop over events
          System.out.printf("\n\nFile %d from %d: wrote %d events....\n\n",nFile,files.length,count);
        }//loop over files
      }//loop over files

      System.out.printf("\n\nWrote %d Events\n\n",count);
      w.close();

    }

    public void skimTwoPionNoElectron(String inDir,String outFile,double beamE, double targetMass, double MMLim, int cutRho, int limEvs){

      File folder = new File(inDir);

      File[] files = folder.listFiles();

      HipoReader rtmp = new HipoReader(inDir+files[0].getName());
      HipoWriter w = HipoWriter.create(outFile, rtmp);

      long count=0, countWrite=0, nFile=0;

      //don't count if set limEvs=-1
      if(limEvs==-1){countWrite=-2;}

      for (File file : files) {
        nFile++;
        if (file.isFile() && countWrite<limEvs) {
          HipoReader r = new HipoReader(inDir+file.getName());
          Event ev = new Event();
          Event ev_out = new Event();

          Bank recpart = r.getBank("REC::Particle");
          Bank reccal = r.getBank("REC::Calorimeter");
          Bank dc = r.getBank("DC::tdc");
          Bank ecal = r.getBank("ECAL::adc");
          Bank ftoftdc = r.getBank("FTOF::tdc");
          Bank ftofadc = r.getBank("FTOF::adc");
          Bank htcc = r.getBank("HTCC::adc");
          Bank rechtcc = r.getBank("REC::Cherenkov");
          Bank recftof = r.getBank("REC::Scintillator");
          Bank rectrack = r.getBank("REC::Track");
          Bank runconfig = r.getBank("RUN::config");

          int nVarsPart=10, nVarsExc=4;
          double[] pim=new double[nVarsPart];
          double[] pip=new double[nVarsPart];
          double[] p=new double[nVarsPart];
          double[] exc=new double[nVarsExc];
          while(r.hasNext() && countWrite<limEvs){

            r.nextEvent(ev);
            ev.read(recpart);
            ev.read(reccal);
            ev.read(rechtcc);
            ev.read(recftof);
            ev.read(rectrack);
            ev.read(runconfig);
            ev.read(ecal);
            ev.read(htcc);
            ev.read(ftoftdc);
            ev.read(ftofadc);
            ev.read(dc);

            int addEv=0;

            for(int row=0;row<recpart.getRows();row++){
              for(int row2=(row+1);row2<recpart.getRows();row2++){
                for(int row3=0;row3<recpart.getRows();row3++){

                  cleanArr(pim,nVarsPart);
                  cleanArr(pip,nVarsPart);
                  cleanArr(p,nVarsPart);
                  cleanArr(exc,nVarsExc);
                  

                  fillRECPart(recpart,rectrack,row,p);
                  fillRECPart(recpart,rectrack,row3,pim);
                  fillRECPart(recpart,rectrack,row2,pip);
                  int hasPCandi=0,hasPim=0,hasPip=0;
                  

                  if(pip[5]==1 && pip[0]==211){
                    hasPip=1;
                  }
                  if(pim[5]==-1 ){ //&& pim[0]==-211
                    hasPim=1;
                  }

                  if(p[5]==1  && p[0]==2212 && row!=row2){ //&& matchel!=-1 el[11] is track chi^2 && el[11]<350 && Math.abs(el[13])<20 && el[12]==6
                    hasPCandi=1;
                  }

                  if(hasPCandi==1 && hasPim==1 && hasPip==1){
                    calcExc(pip,pim, p, exc,beamE,0);
                    //System.out.printf("\nMM %f",exc[1]);
                    if(cutRho==0){
                      if(Math.abs(exc[1])<MMLim){
                        //System.out.printf(" after MM %f \n",exc[1]);
                        addEv=1;
                      }
                    } else{
                      if(Math.abs(exc[1])<MMLim){
                        if(exc[0]>0.6 && exc[0]<0.9){
                          //System.out.printf(" after MM %f \n",exc[1]);
                          addEv=1;
                        }
                      }
                    }//cut on rho or not
                    
                  }//have good candidates

                }//read pi+
              }//read pi-
            }//read el

            if(addEv==1){
              ev_out.reset();
              ev_out.write(recpart);
              ev_out.write(reccal);
              ev_out.write(rechtcc);
              ev_out.write(recftof);
              ev_out.write(rectrack);
              ev_out.write(runconfig);
              ev_out.write(ecal);
              ev_out.write(htcc);
              ev_out.write(ftoftdc);
              ev_out.write(ftofadc);
              ev_out.write(dc);
              w.addEvent(ev_out);
              count++;
              if(limEvs!=-1){
                countWrite++;
              }
            }

          }//loop over events
          System.out.printf("\n\nFile %d from %d: wrote %d events....\n\n",nFile,files.length,count);
        }//loop over files
      }//loop over files

      System.out.printf("\n\nWrote %d Events\n\n",count);
      w.close();

    }

    public void skimTwoPionFTel(String inDir,String outFile,double beamE, double targetMass, double MMUpperLim, int cutRho, int limEvs){

      File folder = new File(inDir);

      File[] files = folder.listFiles();

      HipoReader rtmp = new HipoReader(inDir+files[0].getName());
      HipoWriter w = HipoWriter.create(outFile, rtmp);

      long count=0, countWrite=0, nFile=0;

      //don't count if set limEvs=-1
      if(limEvs==-1){countWrite=-2;}

      for (File file : files) {
        nFile++;
        if (file.isFile() && countWrite<limEvs) {
          HipoReader r = new HipoReader(inDir+file.getName());
          Event ev = new Event();
          Event ev_out = new Event();

          Bank recpart = r.getBank("REC::Particle");
          Bank recftpart = r.getBank("RECFT::Particle");
          Bank reccal = r.getBank("REC::Calorimeter");
          Bank dc = r.getBank("DC::tdc");
          Bank ecal = r.getBank("ECAL::adc");
          Bank ftoftdc = r.getBank("FTOF::tdc");
          Bank ftofadc = r.getBank("FTOF::adc");
          Bank htcc = r.getBank("HTCC::adc");
          Bank rechtcc = r.getBank("REC::Cherenkov");
          Bank recftof = r.getBank("REC::Scintillator");
          Bank rectrack = r.getBank("REC::Track");
          Bank runconfig = r.getBank("RUN::config");

          int nVarsPart=10, nVarsExc=4;
          double[] pim=new double[nVarsPart];
          double[] pip=new double[nVarsPart];
          double[] el=new double[nVarsPart];
          double[] exc=new double[nVarsExc];
          while(r.hasNext() && countWrite<limEvs){

            r.nextEvent(ev);
            ev.read(recpart);
            ev.read(recftpart);
            ev.read(reccal);
            ev.read(rechtcc);
            ev.read(recftof);
            ev.read(rectrack);
            ev.read(runconfig);
            ev.read(ecal);
            ev.read(htcc);
            ev.read(ftoftdc);
            ev.read(ftofadc);
            ev.read(dc);

            int addEv=0;

            if(recftpart.getRows()>0){
              System.out.println("Event with RECFT::Part \n");
              recftpart.show();
              recpart.show();
            }

            for(int row=0;row<recftpart.getRows();row++){
              for(int row2=(row+1);row2<recftpart.getRows();row2++){
                for(int row3=0;row3<recftpart.getRows();row3++){

                  cleanArr(pim,nVarsPart);
                  cleanArr(pip,nVarsPart);
                  cleanArr(el,nVarsPart);
                  cleanArr(exc,nVarsExc);
                  

                  fillRECPart(recftpart,rectrack,row,el);
                  fillRECPart(recftpart,rectrack,row2,pim);
                  fillRECPart(recftpart,rectrack,row3,pip);
                  int hasElCandi=0,hasPim=0,hasPip=0;

                  

                  if(pip[5]==1 && pip[0]==211){
                    hasPip=1;
                  }
                  if(pim[5]==-1 ){ //&& pim[0]==-211
                    hasPim=1;
                  }

                  if(el[5]==-1 && Math.abs(el[4])<2000){ //&& matchel!=-1 el[11] is track chi^2 && el[11]<350 && Math.abs(el[13])<20 && el[12]==6
                    hasElCandi=1;
                  }

                  if(hasElCandi==1 && hasPim==1 && hasPip==1){
                    calcExc(pip,pim, el, exc,beamE,1);
                    //System.out.printf("\nMM %f",exc[1]);
                    if(cutRho==0){
                      if(exc[1]>0.2 && exc[1]<MMUpperLim){
                        //System.out.printf(" after MM %f \n",exc[1]);
                        addEv=1;
                      }
                    } else{
                      if(exc[1]>0.2 && exc[1]<MMUpperLim){
                        if(exc[0]>0.6 && exc[0]<0.9){
                          //System.out.printf(" after MM %f \n",exc[1]);
                          addEv=1;
                        }
                      }
                    }//cut on rho or not
                    
                  }//have good candidates

                }//read pi+
              }//read pi-
            }//read el

            if(addEv==1){
              ev_out.reset();
              ev_out.write(recpart);
              ev_out.write(recftpart);
              ev_out.write(reccal);
              ev_out.write(rechtcc);
              ev_out.write(recftof);
              ev_out.write(rectrack);
              ev_out.write(runconfig);
              ev_out.write(ecal);
              ev_out.write(htcc);
              ev_out.write(ftoftdc);
              ev_out.write(ftofadc);
              ev_out.write(dc);
              w.addEvent(ev_out);
              count++;
              if(limEvs!=-1){
                countWrite++;
              }
            }

          }//loop over events
          System.out.printf("\n\nFile %d from %d: wrote %d events....\n\n",nFile,files.length,count);
        }//loop over files
      }//loop over files

      System.out.printf("\n\nWrote %d Events\n\n",count);
      w.close();

    }

    public void skimKLambda(String inDir,String outFile,double beamE, double targetMass, double MMUpperLim, int cutLambda, int limEvs, int FidStr){

      File folder = new File(inDir);

      File[] files = folder.listFiles();

      HipoReader rtmp = new HipoReader(inDir+files[0].getName());
      HipoWriter w = HipoWriter.create(outFile, rtmp);

      long count=0, countWrite=0, nFile=0;

      //don't count if set limEvs=-1
      if(limEvs==-1){countWrite=-2;}

      for (File file : files) {
        nFile++;
        if (file.isFile() && countWrite<limEvs) {
          HipoReader r = new HipoReader(inDir+file.getName());
          Event ev = new Event();
          Event ev_out = new Event();

          Bank recpart = r.getBank("REC::Particle");
          Bank reccal = r.getBank("REC::Calorimeter");
          Bank dc = r.getBank("DC::tdc");
          Bank ecal = r.getBank("ECAL::adc");
          Bank ftoftdc = r.getBank("FTOF::tdc");
          Bank ftofadc = r.getBank("FTOF::adc");
          Bank htcc = r.getBank("HTCC::adc");
          Bank rechtcc = r.getBank("REC::Cherenkov");
          Bank recftof = r.getBank("REC::Scintillator");
          Bank rectrack = r.getBank("REC::Track");
          Bank runconfig = r.getBank("RUN::config");

          int nVarsPart=10, nVarsExc=4;
          double[] pim=new double[nVarsPart];
          double[] p=new double[nVarsPart];
          double[] el=new double[nVarsPart];
          double[] exc=new double[nVarsExc];
          double[] elEs=new double[9];
          double[] elLs=new double[9];
          double[] pimEs=new double[9];
          double[] pimLs=new double[9];
          double[] pEs=new double[9];
          double[] pLs=new double[9];
          while(r.hasNext() && countWrite<limEvs){

            r.nextEvent(ev);
            ev.read(recpart);
            ev.read(reccal);
            ev.read(rechtcc);
            ev.read(recftof);
            ev.read(rectrack);
            ev.read(runconfig);
            ev.read(ecal);
            ev.read(htcc);
            ev.read(ftoftdc);
            ev.read(ftofadc);
            ev.read(dc);

            int hasTrig=hasTriggerEl(runconfig);

            int addEv=0;

            for(int row=0;row<recpart.getRows();row++){
              for(int row2=(row+1);row2<recpart.getRows();row2++){
                for(int row3=0;row3<recpart.getRows();row3++){

                  cleanArr(pim,nVarsPart);
                  cleanArr(p,nVarsPart);
                  cleanArr(el,nVarsPart);
                  cleanArr(exc,nVarsExc);
                  

                  fillRECPart(recpart,rectrack,row,el);
                  fillRECPart(recpart,rectrack,row2,pim);
                  fillRECPart(recpart,rectrack,row3,p);
                  int hasElCandi=0,hasPim=0,hasP=0;

                  //use el fid to get good sample of negative pions
                  Boolean fid=true;
                  
                  //fids adds extra detector, makes things slow
                  //ignore with fidstrength=-1
                  if(FidStr!=-1){
                    cleanArr(elEs,9);
                    cleanArr(elLs,9);
                    cleanArr(pimEs,9);
                    cleanArr(pimLs,9);
                    cleanArr(pEs,9);
                    cleanArr(pLs,9);
                    getCalInfo(reccal, row, el[1], elEs, elLs);
                    getCalInfo(reccal, row2, pim[1], pimEs, pimLs);
                    getCalInfo(reccal, row3, p[1], pEs, pLs);
                    passFid(elLs,el,FidStr);
                  }
                  

                  if(p[5]==1 && p[0]==2212){
                    hasP=1;
                  }
                  if(pim[5]==-1 ){ //&& pim[0]==-211
                    hasPim=1;
                  }

                  if(el[5]==-1 && fid){ //&& matchel!=-1 el[11] is track chi^2 && el[11]<350 && Math.abs(el[13])<20 && el[12]==6
                    hasElCandi=1;
                  }

                  if(hasElCandi==1 && hasPim==1 && hasP==1 && hasTrig==1){
                    calcExc(p,pim, el, exc,beamE,1);
                    //System.out.printf("\nMM %f",exc[1]);
                    if(cutLambda==0){
                      if(exc[1]>0 && exc[1]<MMUpperLim){
                        //System.out.printf(" after MM %f \n",exc[1]);
                        addEv=1;
                      }
                    } else{
                      if(exc[1]>0 && exc[1]<MMUpperLim){
                        if(exc[0]>0.9 && exc[0]<1.3){
                          //System.out.printf(" after MM %f \n",exc[1]);
                          addEv=1;
                        }
                      }
                    }//cut on rho or not
                    
                  }//have good candidates

                }//read pi+
              }//read pi-
            }//read el

            if(addEv==1){
              ev_out.reset();
              ev_out.write(recpart);
              ev_out.write(reccal);
              ev_out.write(rechtcc);
              ev_out.write(recftof);
              ev_out.write(rectrack);
              ev_out.write(runconfig);
              ev_out.write(ecal);
              ev_out.write(htcc);
              ev_out.write(ftoftdc);
              ev_out.write(ftofadc);
              ev_out.write(dc);
              w.addEvent(ev_out);
              count++;
              if(limEvs!=-1){
                countWrite++;
              }
            }

          }//loop over events
          System.out.printf("\n\nFile %d from %d: wrote %d events....\n\n",nFile,files.length,count);
        }//loop over files
      }//loop over files

      System.out.printf("\n\nWrote %d Events\n\n",count);
      w.close();

    }

    public void skim(String inDir,String outFile){

      File folder = new File(inDir);

      File[] files = folder.listFiles();

      HipoReader rtmp = new HipoReader(inDir+files[0].getName());
      HipoWriter w = HipoWriter.create(outFile, rtmp);

      long count=0;

      for (File file : files) {
        if (file.isFile()) {
          HipoReader r = new HipoReader(inDir+file.getName());
          Event ev = new Event();
          Event ev_out = new Event();

          Bank recpart = r.getBank("REC::Particle");
          Bank reccal = r.getBank("REC::Calorimeter");
          Bank dc = r.getBank("DC::tdc");
          Bank ecal = r.getBank("ECAL::adc");
          Bank ftoftdc = r.getBank("FTOF::tdc");
          Bank ftofadc = r.getBank("FTOF::adc");
          Bank htcc = r.getBank("HTCC::adc");
          Bank rechtcc = r.getBank("REC::Cherenkov");
          Bank recftof = r.getBank("REC::Scintillator");
          Bank rectrack = r.getBank("REC::Track");
          Bank runconfig = r.getBank("RUN::config");

          double[] electron = new double[4];

          double[] photon = new double[4];

          while(r.hasNext()){

            r.nextEvent(ev);
            ev.read(recpart);
            ev.read(reccal);
            ev.read(rechtcc);
            ev.read(recftof);
            ev.read(rectrack);
            ev.read(runconfig);
            ev.read(ecal);
            ev.read(htcc);
            ev.read(ftoftdc);
            ev.read(ftofadc);
            ev.read(dc);

            //works but we also want cases where els are missIDed by recon
            /*
            int elpindex=hasFDPID(recpart,11);
            if(elpindex!=-1){
              if(radiatedPhoton(recpart,elpindex) || hasFDPID(recpart,211)!=-1){
                w.addEvent(ev);
                count++;
              }
            }
            */


            int countNeg=countChargeFD(recpart,-1);
            if(countNeg>0){
              if(hasFDPID(recpart,211)!=-1 || hasPhotonElectronPair(recpart)){
                ev_out.reset();
                ev_out.write(recpart);
                ev_out.write(reccal);
                ev_out.write(rechtcc);
                ev_out.write(recftof);
                ev_out.write(rectrack);
                ev_out.write(runconfig);
                ev_out.write(ecal);
                ev_out.write(htcc);
                ev_out.write(ftoftdc);
                ev_out.write(ftofadc);
                ev_out.write(dc);
                w.addEvent(ev_out);
                count++;
              }
            }
          }//read events
          System.out.printf("\n\nStill reading, so far wrote %d Events....\n\n",count);
        }//check is file
      }//loop over files
      w.close();
      System.out.printf("\n\nWrote %d Events\n\n",count);
    }

    //run with /open src/main/java/j4np/instarec/data/FileSkimmer.java
    //FileSkimmer.main(new String[]{})

    public static void main(String[] args){
        
      System.out.println("\n\n----- FileSkimmer");

      FileSkimmer dp = new FileSkimmer();

      /*String inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      String outName="/w/work/clas12/tyson/data_repo/caos/rga/run_valid_5197.hipo";
      dp.skim(inDir,outName);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005407/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_valid_5407.hipo";
      dp.skim(inDir,outName);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005342/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_valid_5342.hipo";
      dp.skim(inDir,outName);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005418/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_valid_5418.hipo";
      dp.skim(inDir,outName);*/

      //outbending below
      /*inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005442/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_valid_5442.hipo";
      dp.skim(inDir,outName);

      
      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005444/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_valid_5444.hipo";
      dp.skim(inDir,outName);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005543/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_valid_5543.hipo";
      dp.skim(inDir,outName);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005595/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_valid_5595.hipo";
      dp.skim(inDir,outName);*/

      //Twopion

      /*String inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      String outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPion_noPimID_5197.hipo";
      dp.skimTwoPion(inDir,outName,10.6,0.938272,2.1,0,-1,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005407/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPion_noPimID_5407.hipo";
      dp.skimTwoPion(inDir,outName,10.6,0.938272,2.1,0,-1,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005342/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPion_noPimID_5342.hipo";
      dp.skimTwoPion(inDir,outName,10.6,0.938272,2.1,0,-1,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005418/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPion_noPimID_5418.hipo";
      dp.skimTwoPion(inDir,outName,10.6,0.938272,2.1,0,-1,-1);*/

      //outbending below
      /*
      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005442/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPion_5442.hipo";
      dp.skimTwoPion(inDir,outName,10.6,0.938272,2.1,0,-1,-1);
            
      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005444/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPion_5444.hipo";
      dp.skimTwoPion(inDir,outName,10.6,0.938272,2.1,0,-1,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005543/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPion_5543.hipo";
      dp.skimTwoPion(inDir,outName,10.6,0.938272,2.1,0,-1,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005595/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPion_5595.hipo";
      dp.skimTwoPion(inDir,outName,10.6,0.938272,2.1,0,-1,-1);*/

      //twopion, with cuts for negative sample

      //inDir="/volatile/clas12/gavalian/cook/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      //outName="/w/work/clas12/tyson/data_repo/caos/rga/run_trainTwoPion_noPimID_5197.hipo";
      //dp.skimTwoPion(inDir,outName,10.6,0.938272,1.1,0,-1,-1);

      //KLambda

      /*String inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      String outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validKLambda_5197.hipo";
      dp.skimKLambda(inDir,outName,10.6,0.938272,2.1,0,-1,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005407/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validKLambda_5407.hipo";
      dp.skimKLambda(inDir,outName,10.6,0.938272,2.1,0,-1,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005342/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validKLambda_5342.hipo";
      dp.skimKLambda(inDir,outName,10.6,0.938272,2.1,0,-1,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005418/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validKLambda_5418.hipo";
      dp.skimKLambda(inDir,outName,10.6,0.938272,2.1,0,-1,-1);*/

      //large
      /*String inDir="/volatile/clas12/gavalian/cook/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      String outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validKLambda_5197_large.hipo";
      dp.skimKLambda(inDir,outName,10.6,0.938272,1.0,0,-1,-1);*/


      //TwoPionFTel

      /*String inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      String outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPionFTel_5197.hipo";
      dp.skimTwoPionFTel(inDir,outName,10.6,0.938272,2.1,0,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005407/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPionFTel_5407.hipo";
      dp.skimTwoPionFTel(inDir,outName,10.6,0.938272,2.1,0,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005342/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPionFTel_5342.hipo";
      dp.skimTwoPionFTel(inDir,outName,10.6,0.938272,2.1,0,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005418/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPionFTel_5418.hipo";
      dp.skimTwoPionFTel(inDir,outName,10.6,0.938272,2.1,0,-1);*/

      //large
      //String inDir="/volatile/clas12/gavalian/cook/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      //String outName="/w/work/clas12/tyson/data_repo/caos/rga/run_trainTwoPionFTel_5197_large.hipo";
      //dp.skimTwoPionFTel(inDir,outName,10.6,0.938272,1.0,0,-1);


      String inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      String outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPionNoEl_5197.hipo";
      dp.skimTwoPionNoElectron(inDir,outName,10.6,0.938272,0.5,0,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005407/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPionNoEl_5407.hipo";
      dp.skimTwoPionNoElectron(inDir,outName,10.6,0.938272,0.5,0,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005342/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPionNoEl_5342.hipo";
      dp.skimTwoPionNoElectron(inDir,outName,10.6,0.938272,0.5,0,-1);

      inDir="/volatile/clas12/users/caot/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005418/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_validTwoPionNoEl_5418.hipo";
      dp.skimTwoPionNoElectron(inDir,outName,10.6,0.938272,0.5,0,-1);

      //large
      inDir="/volatile/clas12/gavalian/cook/experiment/aiAssistedPlusTracking_iss471/rga_fall2018/full/recon/005197/"; // used to change output path of plots (eg adding _NoFiducialCuts)
      outName="/w/work/clas12/tyson/data_repo/caos/rga/run_trainTwoPionNoEl_5197.hipo";
      dp.skimTwoPionNoElectron(inDir,outName,10.6,0.938272,0.5,0,-1);
        
    }

}