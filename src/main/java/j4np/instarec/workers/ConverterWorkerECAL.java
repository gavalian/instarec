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

/**
 *
 * @author tyson
 */
public class ConverterWorkerECAL extends DataWorker {
    Schema ecalSchema = null;
    public ConverterWorkerECAL(){
        //DC::tdc/20600/12,"sector/B,layer/B,component/S,order/B,TDC/I";
        SchemaBuilder b = new SchemaBuilder("ECAL::adc",20700,11);
        ecalSchema = b.addEntry("sector", "B", "")
                .addEntry("layer", "B", "")
                .addEntry("component", "S", "")
                .addEntry("order", "B", "")
                .addEntry("ADC", "I", "")
                .addEntry("time", "F", "")
                .addEntry("ped", "S", "").build();
        
        //ecalSchema.show();
    }


    @Override
    public boolean init(DataSource src) {
        //throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
        return true;
    }

    @Override
    public void execute(DataEvent event) {
        Bank b = new Bank(ecalSchema,1024);
        ((Event) event).read(b);
        Leaf adc = new Leaf(42,12,"bbsf",4096);

        adc.setRows(b.getRows());
        //System.out.printf(" rows = %d , leaf = %d\n",b.getRows(), adc.getRows());
        for(int r = 0; r < b.getRows(); r++){
          //System.out.printf("sector %d layer %d comp %d adc %d \n",b.getByte(0, r),b.getByte(1, r),b.getShort(2, r),b.getInt(4, r));
            adc.putByte(0, r, b.getByte(0, r));
            adc.putByte(1, r, b.getByte(1, r));
            adc.putShort(2, r, b.getShort(2, r));
            adc.putFloat(3,  r, b.getInt(4, r));
        }
        ((Event) event).write(adc);
    }
    
}
