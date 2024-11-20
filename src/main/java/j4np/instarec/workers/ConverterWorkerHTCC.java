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
import j4np.hipo5.io.HipoReader;

/**
 *
 * @author tyson
 */
public class ConverterWorkerHTCC extends DataWorker {
    Schema htccSchema = null;
    public ConverterWorkerHTCC(){
        //DC::tdc/20600/12,"sector/B,layer/B,component/S,order/B,TDC/I";
        SchemaBuilder b = new SchemaBuilder("HTCC::adc",21500,11);
        htccSchema = b.addEntry("sector", "B", "")
                .addEntry("layer", "B", "")
                .addEntry("component", "S", "")
                .addEntry("order", "B", "")
                .addEntry("ADC", "I", "")
                .addEntry("time", "F", "")
                .addEntry("ped", "S", "").build();
        
        //dcSchema.show();
    }


    @Override
    public boolean init(DataSource src) {
        //throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
        return true;
    }

    @Override
    public void execute(DataEvent event) {
        Bank b = new Bank(htccSchema,4096);
        ((Event) event).read(b);
        Leaf adc = new Leaf(42,13,"i8f",4096);
        
        //intialise bank, will be easier to always have all 6 sectors
        adc.setRows(6);
        for(int i=0;i<6;i++){
          adc.putInt(0, i, i+1);
          for (int j=1;j<9;j++){
            adc.putFloat(j, i, 0);
          }
        }

        for(int r = 0; r < b.getRows(); r++){
          if(b.getInt(4, r)>0.){
            int index = ((b.getByte(1, r) - 1) * 4 + b.getShort(2, r)); //-1 to start from 0-7
            adc.putFloat(index, b.getByte(0, r)-1, b.getInt(4, r));
          }
        }

        ((Event) event).write(adc);
    }
    
}
