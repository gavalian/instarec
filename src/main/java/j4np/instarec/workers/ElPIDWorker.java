/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.workers;

import j4np.data.base.DataActor;
import j4np.data.base.DataEvent;
import j4np.data.base.DataSource;
import j4np.data.base.DataWorker;
import j4np.hipo5.data.Bank;
import j4np.hipo5.data.Event;
import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Schema;
import j4np.hipo5.data.Schema.SchemaBuilder;
import j4np.instarec.utils.NeuralModel;

/**
 *
 * @author tyson
 */
public class ElPIDWorker extends DataWorker {
    int nPredictedOutput=2;
    double threshold_on_response=0.1;
    NeuralModel[] ElPIDModels = new NeuralModel[6];

    //to fill array used for prediction
    //and output array, need to know
    //where to start putting variables
    //call this offsets
    int stripsoffset=3,stripsoffsetout=18;
    int Dsoffset=12, htccoffset=27, htccoffsetout=30;
    int wiresoffset=21,wiresoffsetout=12,wiresoffsetin=17;
    int ftofoffsetout=27;

    public ElPIDWorker(String networkPath,double th){
      threshold_on_response=th;
      for(int i=1;i<7;i++){
        ElPIDModels[i-1]=NeuralModel.jsonFile(networkPath+"_sector"+String.valueOf(i)+".json");
      }
    }


    @Override
    public boolean init(DataSource src) {
        //throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
        return true;
    }

    public void readTrackBank(Leaf trackbank, Leaf partout, float[] nInVars, short row){
      //add ps and vs from ai tracking in output
      for (int j=6;j<12;j++){ partout.putFloat(j,row,(float)trackbank.getDouble(j-1,row));}
      //fill wires, take average between two slopes
      for (int j=0;j<6;j++){
        nInVars[j+wiresoffset]=((float)trackbank.getDouble(j+wiresoffsetin,row)+(float)trackbank.getDouble(j+wiresoffsetin+6,row))/2;
        partout.putFloat(j+wiresoffsetout,row, nInVars[j+wiresoffset]);
      }
    }

    public double readECALClusters(Leaf ECALclusters, Leaf partout, float[] nInVars, short row){
      double sumE=0;
      //init sums of edep to 0
      for(int j=0;j<3;j++){nInVars[j]=0;}
      for(int j=0;j<ECALclusters.getRows();j++){
        //match cluster to track
        if(ECALclusters.getShort(0,j)==row){
          if(ECALclusters.getInt(2,j)!=0){

            nInVars[ECALclusters.getInt(2,j)-1+stripsoffset]=(float)ECALclusters.getDouble(3,j);
            nInVars[ECALclusters.getInt(2,j)-1+Dsoffset]=ECALclusters.getInt(5,j);
            partout.putFloat(ECALclusters.getInt(2,j)-1+stripsoffsetout,row,(float)ECALclusters.getDouble(3,j));
            sumE+=ECALclusters.getDouble(4,j);

            //sum energy in PCAL, ECIN, ECOUT
            if(ECALclusters.getInt(2,j)<4){
              nInVars[0]+=ECALclusters.getDouble(4,j);
            } else if(ECALclusters.getInt(2,j)>=4 && ECALclusters.getInt(2,j)<7){
              nInVars[1]+=ECALclusters.getDouble(4,j);
            } else if(ECALclusters.getInt(2,j)>=7){
              nInVars[2]+=ECALclusters.getDouble(4,j);
            }
          }
        }
      }
      return sumE;
    }

    public float readHTCCadc(Leaf adc, Leaf partout, float[] nInVars, short row, short sector){

      float sumHTCC=0, sumHTCCBefore=0,sumHTCCAfter=0;

      int sectorBefore=sector-1;
      int sectorAfter=sector+1;
      if(sector==1){sectorBefore=6;}
      else if(sector==6){sectorAfter=1;}

      for(int k = 0; k < adc.getRows(); k++){
        int detector_type = adc.getInt(0, k); //no get byte function??
        int sect = adc.getInt(1, k); //no get byte function??
        int layer = adc.getInt(2, k);
        int strip = adc.getInt(3, k);
        int ADC = adc.getInt(5, k);
        if(ADC>0. && detector_type==15){
          int index = ((layer - 1) * 4 + strip) - 1; //-1 for 0-7, not for 1-8

          if(sect==sector){
            nInVars[index+htccoffset]=ADC;
            sumHTCC+=ADC;
          } else if(sect==sectorBefore){
            nInVars[index+htccoffset+8]=ADC;
            sumHTCCBefore+=ADC;
          } else if(sect==sectorAfter){
            nInVars[index+htccoffset+16]=ADC;
            sumHTCCAfter+=ADC;
          }
        }
      }

      partout.putFloat(htccoffsetout,row,sumHTCC);
      partout.putFloat(htccoffsetout+1,row,sumHTCCBefore);
      partout.putFloat(htccoffsetout+2,row,sumHTCCAfter);

      return (sumHTCC+sumHTCCAfter+sumHTCCBefore);
    }

    public float readFTOFClusters(Leaf FTOFclusters, Leaf partout, float[] nInVars, short row){
      float FTOFpredPath=0;
      //fill FTOF info
      for(int j=0;j<FTOFclusters.getRows();j++){
        if(FTOFclusters.getShort(0,j)==row && FTOFclusters.getInt(2,j)==2){
          FTOFpredPath=(float)FTOFclusters.getDouble(4,j);
          partout.putFloat(ftofoffsetout,row,FTOFpredPath); //pred path
          partout.putFloat(ftofoffsetout+1,row,(float)FTOFclusters.getDouble(6,j)); //pred time
          partout.putFloat(ftofoffsetout+2,row,(float)FTOFclusters.getDouble(3,j)); //pred comp
        }
      }
      return FTOFpredPath;
    }

    public void predictPID(Leaf partout, float[] nInVars, float[] pred_pid,int charge, int sector, short row, double sumE, float sumHTCC, float FTOFpredPath){
      if(charge==-1 && sector>0){
        ElPIDModels[sector-1].predict(nInVars, pred_pid);
        partout.putFloat(2,row,pred_pid[1]);
        //for e- PID require good pred
        //but also some sanity checks:
        //non empty HTCC, ECAL, and track hit FTOF (eg good track)
        if(pred_pid[1]>threshold_on_response && FTOFpredPath>0 && sumHTCC>0  && sumE>0){
          partout.putInt(1,row,11);
        } else {
          //dummy PID, need worker for non el PID
          partout.putInt(1,row,-211); 
        }

      } else{
        //dummy PID, need worker for non el PID
        partout.putInt(1,row,211); 
        partout.putFloat(2,row,(float)0.5);
      }
    }
    
 
    @Override
    public void execute(DataEvent event) {
        Leaf trackbank = new Leaf(32000,1,"i",4096);
        ((Event) event).read(trackbank);
        Leaf ECALclusters = new Leaf(32200,2,"i",4096);
        ((Event) event).read(ECALclusters);
        Leaf adc = new Leaf(42,12,"i",4096); //read htcc adc from here
        ((Event) event).read(adc);
        Leaf FTOFclusters = new Leaf(32200,3,"i",4096);
        ((Event) event).read(FTOFclusters);
        //pindex, pid, pid prob, sector, charge, beta, pxpypz, vxvyvz
        //6xwires, 9xec clusters, ftof (layer 2) path/time/component
        //HTCC sum ADC (same sector, before, after)
        Leaf partout = new Leaf(32200,1,"sifssf3f3f6f9f3ffff",4096);

        partout.setRows(trackbank.getRows()); 
        float[] nInVars=new float[51];
        float[] predicted_pid=new float[2];

        for(short i=0;i<trackbank.getRows();i++){
          short charge=trackbank.getShort(3,i);
          short sector=trackbank.getShort(2,i);

          //redundant but linking btw track and particle might become important
          partout.putShort(0,i,i); 
          
          partout.putFloat(5,i,1); //dummy beta, need other worker
          partout.putShort(3,i,sector);
          partout.putShort(4,i,charge);
          
          readTrackBank(trackbank,partout,nInVars,i);

          double sumE=readECALClusters(ECALclusters,partout,nInVars,i);

          float FTOFpredPath=readFTOFClusters(FTOFclusters,partout,nInVars,i);
          
          float sumHTCC=readHTCCadc(adc,partout,nInVars,i, sector);

          predictPID(partout, nInVars, predicted_pid, charge, sector, i, sumE, sumHTCC, FTOFpredPath);

        }

        ((Event) event).write(partout);
    }
    
}
