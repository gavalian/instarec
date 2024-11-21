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
public class ConverterWorkerParticleCFTraining extends DataWorker {
    Schema recparticleSchema = null;
    Schema reccalorimeterSchema = null;
    Schema recftofSchema = null;
    Schema rechtccSchema = null;
    double lim_px=0.2,lim_py=0.2,lim_pz=0.2;

    public ConverterWorkerParticleCFTraining(){
        //DC::tdc/20600/12,"sector/B,layer/B,component/S,order/B,TDC/I";
        SchemaBuilder brp = new SchemaBuilder("REC::Particle",300,31);
        recparticleSchema = brp.addEntry("pid", "I","")
                      .addEntry("px", "F","")
                      .addEntry("py", "F","")
                      .addEntry("pz", "F","")
                      .addEntry("vx", "F","")
                      .addEntry("vy", "F","")
                      .addEntry("vz", "F","")
                      .addEntry("vt", "F","")
                      .addEntry("charge", "B","")
                      .addEntry("beta", "F","")
                      .addEntry("chi2pid", "F","")
                      .addEntry("status", "S","").build();

        SchemaBuilder brc = new SchemaBuilder("REC::Calorimeter",300,32);
        reccalorimeterSchema = brc.addEntry("index", "S","")
                      .addEntry("pindex", "S","")
                      .addEntry("detector", "B","")
                      .addEntry("sector", "B","")
                      .addEntry("layer", "B","")
                      .addEntry("energy", "F","")
                      .addEntry("time", "F","") 
                      .addEntry("path", "F","")
                      .addEntry("chi2", "F","")
                      .addEntry("x", "F","")
                      .addEntry("y", "F","")
                      .addEntry("z", "F","")
                      .addEntry("hx", "F","")
                      .addEntry("hy", "F","")
                      .addEntry("hz", "F","")
                      .addEntry("lu", "F","")
                      .addEntry("lv", "F","")
                      .addEntry("lw", "F","")
                      .addEntry("du", "F","")
                      .addEntry("dv", "F","")
                      .addEntry("dw", "F","")
                      .addEntry("m2u", "F","")
                      .addEntry("m2v", "F","")
                      .addEntry("m2w", "F","")
                      .addEntry("m3u", "F","")
                      .addEntry("m3v", "F","")
                      .addEntry("m3w", "F","")
                      .addEntry("status", "S","").build();

        SchemaBuilder bftof = new SchemaBuilder("REC::Scintillator",300,35);
        recftofSchema = bftof.addEntry("index", "S","")
                      .addEntry("pindex","S","")
                      .addEntry("detector","B","")
                      .addEntry("sector","B","")
                      .addEntry("layer", "B","")
                      .addEntry("component","S","")
                      .addEntry("energy","F","")
                      .addEntry("time","F","")
                      .addEntry("path","F","")
                      .addEntry("chi2","F","")
                      .addEntry("x","F","")
                      .addEntry("y","F","")
                      .addEntry("z","F","")
                      .addEntry("hx","F","")
                      .addEntry("hy","F","")
                      .addEntry("hz","F","")
                      .addEntry("status","S","").build();

        SchemaBuilder bhtcc = new SchemaBuilder("REC::Cherenkov",300,33);
        rechtccSchema = bhtcc.addEntry("index", "S","")
                      .addEntry("pindex","S","")
                      .addEntry("detector","B","")
                      .addEntry("sector","B","")
                      .addEntry("nphe","F","")
                      .addEntry("time","F","")
                      .addEntry("path","F","")
                      .addEntry("chi2","F","")
                      .addEntry("x","F","")
                      .addEntry("y","F","")
                      .addEntry("z","F","")
                      .addEntry("dtheta","F","")
                      .addEntry("dphi","F","")
                      .addEntry("status","S","").build();
        
    }

    public int fillRECParticleInfo(Leaf partout,Bank PartBank, int row, int pred_charge, double pred_px, double pred_py, double pred_pz) {
        
      int best_pid=1;
      short best_pindex=999;
      float best_resz=999,best_resx=999,best_resy=999, best_beta=999;
      for (short i = 0; i < PartBank.getRows(); i++) {
        int pid = PartBank.getInt("pid", i);
        int status = PartBank.getShort("status", i);
        int charge = PartBank.getByte("charge", i);
        float pz = PartBank.getFloat("pz", i);
        float px = PartBank.getFloat("px", i);
        float py = PartBank.getFloat("py", i);
        float vz = PartBank.getFloat("vz", i);
        float chi2pid = PartBank.getFloat("chi2pid", i);
        float beta = PartBank.getFloat("beta", i);
        //only interested in FD particles for now
        if (Math.abs(status) >= 2000 && Math.abs(status) < 4000) {
          //want  good particles, only using this to create training sample
          if(vz<12 && vz>-13 && Math.abs(chi2pid)<5){
            float resz=(float)pred_pz-pz;
            float resx=(float)pred_px-px;
            float resy=(float)pred_py-py;
            if( Math.abs(resz)<best_resz && resz<lim_pz && Math.abs(resx)<best_resx && resx<lim_px && Math.abs(resy)<best_resy && resy<lim_py&& pred_charge==charge){ 
              best_pindex=i;
              best_pid=pid;
              best_beta=beta;
              best_resx=resx;
              best_resy=resy;
              best_resz=resz;
            }
          }
        }
      }

      partout.putShort(0,row,best_pindex);
      partout.putInt(1,row,best_pid);
      partout.putFloat(2,row,1);
      partout.putFloat(5,row,best_beta);

      return best_pindex;
    }

    
    public void fillRECCalorimeterInfo(Leaf partout,Bank ECAL_Bank,int row, int pindex,short pred_sector) {
      
      for (int k = 0; k < ECAL_Bank.getRows(); k++) {
        short i = ECAL_Bank.getShort("pindex", k);
        byte sector = ECAL_Bank.getByte("sector", k);
        float lu=ECAL_Bank.getFloat("lu",k);
        float lv=ECAL_Bank.getFloat("lv",k);
        float lw=ECAL_Bank.getFloat("lw",k);
        byte layer=ECAL_Bank.getByte("layer",k);
        if (i == pindex && pred_sector==sector) {
          //Cal_index=index;
          if(layer==1){
            partout.putFloat(18,row,lu);
            partout.putFloat(19,row,lv);
            partout.putFloat(20,row,lw);
          } else if(layer==4){
            partout.putFloat(21,row,lu);
            partout.putFloat(22,row,lv);
            partout.putFloat(23,row,lw);
          } else if(layer==7){
            partout.putFloat(24,row,lu);
            partout.putFloat(25,row,lv);
            partout.putFloat(26,row,lw); 
          }
        }
      }
      
    }

    public void fillRECHTCCInfo(Leaf partout,Bank HTCCBank,int row, int pindex) {
      for (int k = 0; k < HTCCBank.getRows(); k++) {
        short i = HTCCBank.getShort("pindex", k);
        float nphe = HTCCBank.getFloat("nphe", k);
        if(i==pindex){
          partout.putFloat(30,row,nphe);
        }
      }
    }

    public void fillRECFTOFInfo(Leaf partout,Bank Scint_Bank,int row, int pindex,short pred_sector){
      //Scint_Bank.show();
      for (int k = 0; k < Scint_Bank.getRows(); k++) {
          short i = Scint_Bank.getShort("pindex", k);
          byte sector = Scint_Bank.getByte("sector", k);
          byte layer=Scint_Bank.getByte("layer",k);
          short component=Scint_Bank.getShort("component",k);
          float path = Scint_Bank.getFloat("path", k);
          float time = Scint_Bank.getFloat("time", k);
          if (i == pindex && sector==pred_sector) {
            if(layer==2){
              partout.putFloat(27,row,path);
              partout.putFloat(28,row,time);
              partout.putFloat(29,row,component);
            } 
          }
      }
    }


    @Override
    public boolean init(DataSource src) {
        //throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
        return true;
    }

    @Override
    public void execute(DataEvent event) {
        Bank bpart = new Bank(recparticleSchema,1024);
        Bank bcal = new Bank(reccalorimeterSchema,1024);
        Bank bftof = new Bank(recftofSchema,1024);
        Bank bhtcc = new Bank(rechtccSchema,1024);
        ((Event) event).read(bpart);
        ((Event) event).read(bcal);
        ((Event) event).read(bftof);
        ((Event) event).read(bhtcc);
        Leaf track = new Leaf(32000,1,"i",4096);
        ((Event) event).read(track);
        //pindex, pid, pid prob, sector, charge, beta, pxpypz, vxvyvz
        //6xwires, 9xec clusters, ftof (layer 2) path/time/component, HTCC nphe
        Leaf partout = new Leaf(42,15,"sifssf3f3f6f9f3ff",4096);

        /*bpart.show();
        bcal.show();
        bftof.show();
        bhtcc.show();*/

        partout.setRows(track.getRows());
        for(int i=0;i<track.getRows();i++){

          partout.putShort(3,i,track.getShort(2,i));
          partout.putShort(4,i,track.getShort(3,i));
          //add ps and vs from ai tracking
          for (int j=6;j<12;j++){ partout.putFloat(j,i,(float)track.getDouble(j-1,i));}
          //fill wires, take average between two slopes
          for (int j=12;j<18;j++){partout.putFloat(j,i, ((float)track.getDouble(j+5,i)+(float)track.getDouble(j+11,i))/2);}

          int rec_pindex=fillRECParticleInfo(partout, bpart,i,track.getShort(3,i),track.getDouble(5,i),track.getDouble(6,i),track.getDouble(7,i));
          fillRECCalorimeterInfo(partout,bcal,i,rec_pindex,track.getShort(2,i));
          fillRECFTOFInfo(partout,bftof,i,rec_pindex,track.getShort(2,i));
          fillRECHTCCInfo(partout,bhtcc,i,rec_pindex);

        }
        
       
        ((Event) event).write(partout);
    }
    
}
