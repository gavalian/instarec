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
public class ClusterFinderWorkerFTOF extends DataWorker {
    Schema ftofSchema = null;
    int nPredictedOutput=11;
    int clsize=1;
    NeuralModel[] cfModels_negatives = new NeuralModel[6];
    NeuralModel[] cfModels_positives = new NeuralModel[6];
    public ClusterFinderWorkerFTOF(String networkPath){
        //DC::tdc/20600/12,"sector/B,layer/B,component/S,order/B,TDC/I";
        SchemaBuilder b = new SchemaBuilder("FTOF::adc",20700,11);
        ftofSchema = b.addEntry("sector", "B", "")
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
        
        //ftofSchema.show();
    }


    @Override
    public boolean init(DataSource src) {
        //throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
        return true;
    }

    public void fillClusterBank(Leaf clusters,Bank FTOF_Bank,float[] pred_cf,short pred_sector, short row){

      float[] timecs= new float[2];
      float[] energycs= new float[2];

      for(int k = 0; k < FTOF_Bank.getRows(); k++){
            
        byte sect = FTOF_Bank.getByte("sector", k);
        byte layer = FTOF_Bank.getByte("layer", k);
        byte order = FTOF_Bank.getByte("order", k);
        short strip = FTOF_Bank.getShort("component", k);
        int ADC = FTOF_Bank.getInt("ADC", k);
        float time = FTOF_Bank.getFloat("time", k);

        if(time>0.0 && sect==pred_sector && layer==2 && order<2){ //orders of 10, 20 ? what does that mean??
          //allow some tolerance for the prediction to be wrong
          if (strip > (Math.round(pred_cf[9]) - clsize) && strip < (Math.round(pred_cf[9]) + clsize)) {
            timecs[order]=time;
            energycs[order]=ADC;
          }
        }

      }

      clusters.putShort(0, row, row);
      clusters.putShort(1, row, pred_sector);
      clusters.putInt(2, row, 2); //only using layer 2 atm
      clusters.putFloat(3,row,pred_cf[9]);
      clusters.putFloat(4,row,pred_cf[10]);
      clusters.putFloat(5,row,(energycs[0]+energycs[1])/2);
      clusters.putFloat(6,row,(timecs[0]+timecs[1])/2);

    }

    @Override
    public void execute(DataEvent event) {
        Bank b = new Bank(ftofSchema,1024);
        ((Event) event).read(b);
        Leaf trackbank = new Leaf(32000,1,"i",4096);
        ((Event) event).read(trackbank);
        //track bank row, sector, layer,  component, path, av energy btw two layers, av time btw two layers
        Leaf clusters = new Leaf(42,14,"ssiffff",4096);

        float[] track= new float[6];
        float[] predicted_clusterpos = new float[nPredictedOutput];
        clusters.setRows(trackbank.getRows()); //1 layer per track, could expand to other ftof layers
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
