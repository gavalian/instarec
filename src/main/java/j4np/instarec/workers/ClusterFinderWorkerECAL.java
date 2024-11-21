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
public class ClusterFinderWorkerECAL extends DataWorker {
    Schema ecalSchema = null;
    int nPredictedOutput=11;
    int clsize=4;
    NeuralModel[] cfModels_negatives = new NeuralModel[6];
    NeuralModel[] cfModels_positives = new NeuralModel[6];
    public ClusterFinderWorkerECAL(String networkPath){
        //DC::tdc/20600/12,"sector/B,layer/B,component/S,order/B,TDC/I";
        SchemaBuilder b = new SchemaBuilder("ECAL::adc",20700,11);
        ecalSchema = b.addEntry("sector", "B", "")
                .addEntry("layer", "B", "")
                .addEntry("component", "S", "")
                .addEntry("order", "B", "")
                .addEntry("ADC", "I", "")
                .addEntry("time", "F", "")
                .addEntry("ped", "S", "").build();

        for(int i=1;i<7;i++){
          cfModels_negatives[i-1]=NeuralModel.jsonFile(networkPath+"_negatives_sector"+String.valueOf(i)+".json");
          cfModels_positives[i-1]=NeuralModel.jsonFile(networkPath+"_positives_sector"+String.valueOf(i)+".json");
        }

        //ecalSchema.show();
    }


    @Override
    public boolean init(DataSource src) {
        //throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
        return true;
    }

    public void fillClusterBank(Leaf clusters,Bank ECAL_Bank,float[] pred_cf,short pred_sector, short row){

      int[] DU= new int[9];
      float[] energy= new float[9];
      for (int i=0;i<9;i++){
        DU[i]=0;
        energy[i]=0;
      }

      for(int k = 0; k < ECAL_Bank.getRows(); k++){
            
        byte sect = ECAL_Bank.getByte("sector", k);
        byte layer = ECAL_Bank.getByte("layer", k);
        short strip = ECAL_Bank.getShort("component", k);
        int ADC = ECAL_Bank.getInt("ADC", k);

        if(ADC>0.0 && sect==pred_sector){
          if (strip > (Math.round(pred_cf[layer-1]) - clsize) && strip < (Math.round(pred_cf[layer-1]) + clsize)) {
            DU[layer-1]++;
            energy[layer-1]+= ADC;
          }
        }
      }

      for(short lay=0;lay<9;lay++){
        clusters.putShort(0, row+lay, row);
        clusters.putShort(1, row+lay, pred_sector);
        clusters.putInt(2, row+lay, lay+1);
        clusters.putFloat(3,row+lay,pred_cf[lay]);
        clusters.putFloat(4,row+lay,energy[lay]);
        clusters.putInt(5,row+lay,DU[lay]);
      }
    }
 
    @Override
    public void execute(DataEvent event) {
        Bank b = new Bank(ecalSchema,1024);
        ((Event) event).read(b);
        Leaf trackbank = new Leaf(32000,1,"i",4096);
        ((Event) event).read(trackbank);
        //track bank row, sector, ecal layer, strip, energy in cluster, n strips fire in cluster
        Leaf clusters = new Leaf(42,12,"ssiffi",4096);

        float[] track= new float[6];
        float[] predicted_clusterpos = new float[nPredictedOutput];
        clusters.setRows(trackbank.getRows()*9); //9 layers per track
        for(short i=0;i<trackbank.getRows();i++){
          //get track
          //take average between prediction with both slopes
          for (int j=0;j<6;j++){track[j]=((float)trackbank.getDouble(j+17,i)+(float)trackbank.getDouble(j+23,i))/2;}

          //find cluster positions
          if(trackbank.getShort(3,i)==-1){
            cfModels_negatives[trackbank.getShort(2,i)-1].predict(track, predicted_clusterpos);
          } else if(trackbank.getShort(3,i)==1){
            cfModels_positives[trackbank.getShort(2,i)-1].predict(track, predicted_clusterpos);
          }

          fillClusterBank(clusters,b,predicted_clusterpos,trackbank.getShort(2,i),i);
          

        }

        ((Event) event).write(clusters);
    }
    
}
