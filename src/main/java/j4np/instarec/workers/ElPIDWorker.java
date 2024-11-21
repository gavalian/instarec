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

    
 
    @Override
    public void execute(DataEvent event) {
        Leaf trackbank = new Leaf(32000,1,"i",4096);
        ((Event) event).read(trackbank);
        Leaf ecalbank = new Leaf(42,12,"i",4096);
        ((Event) event).read(ecalbank);
        Leaf htccbank = new Leaf(42,13,"i",4096);
        ((Event) event).read(htccbank);
        Leaf ftofbank = new Leaf(42,14,"i",4096);
        ((Event) event).read(ftofbank);
        //pindex, pid, pid prob, sector, charge, beta, pxpypz, vxvyvz
        //6xwires, 9xec clusters, ftof (layer 2) path/time/component, HTCC sum ADC
        Leaf partout = new Leaf(32000,2,"sifssf3f3f6f9f3ff",4096);

        partout.setRows(trackbank.getRows()); 
        float[] nInVars=new float[35];
        float[] predicted_pid=new float[2];

        //to fill array used for prediction
        //and output array, need to know
        //where to start putting variables
        //call this offsets
        int stripsoffset=3,stripsoffsetout=18;
        int Dsoffset=12, htccoffset=27, htccoffsetout=30;
        int wiresoffset=21,wiresoffsetout=12,wiresoffsetin=17;
        int ftofoffsetout=27;

        for(short i=0;i<trackbank.getRows();i++){
          short charge=trackbank.getShort(3,i);
          short sector=trackbank.getShort(2,i);

          //redundant but linking btw track and particle might become important
          partout.putShort(0,i,i); 
          
          partout.putFloat(5,i,1); //dummy beta, need other worker
          partout.putShort(3,i,sector);
          partout.putShort(4,i,charge);
          //add ps and vs from ai tracking in output
          for (int j=6;j<12;j++){ partout.putFloat(j,i,(float)trackbank.getDouble(j-1,i));}

          double sumE=0;
          for(int j=0;j<3;j++){nInVars[j]=0;}
          for(int j=0;j<ecalbank.getRows();j++){
            //match cluster to track
            if(ecalbank.getShort(0,j)==i){
              if(ecalbank.getInt(2,j)!=0){
                nInVars[ecalbank.getInt(2,j)-1+stripsoffset]=(float)ecalbank.getDouble(3,j);
                nInVars[ecalbank.getInt(2,j)-1+Dsoffset]=ecalbank.getInt(5,j);
                partout.putFloat(ecalbank.getInt(2,j)-1+stripsoffsetout,i,(float)ecalbank.getDouble(3,j));
                sumE+=ecalbank.getDouble(4,j);
                //sum energy in PCAL, ECIN, ECOUT
                if(ecalbank.getInt(2,j)<4){
                  nInVars[0]+=ecalbank.getDouble(4,j);
                } else if(ecalbank.getInt(2,j)>=4 && ecalbank.getInt(2,j)<7){
                  nInVars[1]+=ecalbank.getDouble(4,j);
                } else if(ecalbank.getInt(2,j)>=7){
                  nInVars[2]+=ecalbank.getDouble(4,j);
                }
              }
            }
          }

          float FTOFpredPath=0;
          //fill FTOF info
          for(int j=0;j<ftofbank.getRows();j++){
            if(ftofbank.getShort(0,j)==i && ftofbank.getInt(2,j)==2){
              FTOFpredPath=(float)ftofbank.getDouble(4,j);
              partout.putFloat(ftofoffsetout,i,FTOFpredPath); //pred path
              partout.putFloat(ftofoffsetout+1,i,(float)ftofbank.getDouble(6,j)); //pred time
              partout.putFloat(ftofoffsetout+2,i,(float)ftofbank.getDouble(3,j)); //pred comp
            }
          }
          
          //fill wires, take average between two slopes
          for (int j=0;j<6;j++){
            nInVars[j+wiresoffset]=((float)trackbank.getDouble(j+wiresoffsetin,i)+(float)trackbank.getDouble(j+wiresoffsetin+6,i))/2;
            partout.putFloat(j+wiresoffsetout,i, nInVars[j+wiresoffset]);
          }

          float sumHTCC=0;
          for(int j=0;j<8;j++){
            nInVars[j+htccoffset]=(float)htccbank.getDouble(j+1,sector-1);
            sumHTCC+=htccbank.getDouble(j+1,sector-1);
          }

          partout.putFloat(htccoffsetout,i,sumHTCC);

          if(charge==-1 && sector>0){
            ElPIDModels[sector-1].predict(nInVars, predicted_pid);
            partout.putFloat(2,i,predicted_pid[1]);
            //for e- PID require good pred
            //but also some sanity checks:
            //non empty HTCC, ECAL, and track hit FTOF (eg good track)
            if(predicted_pid[1]>threshold_on_response && FTOFpredPath>0 && sumHTCC>0  && sumE>0){
              partout.putInt(1,i,11);
            } else {
              //dummy PID, need worker for non el PID
              partout.putInt(1,i,211); 
            }

          } else{
            //dummy PID, need worker for non el PID
            partout.putInt(1,i,211); 
            partout.putFloat(2,i,(float)0.5);
          }
          

        }

        ((Event) event).write(partout);
    }
    
}
