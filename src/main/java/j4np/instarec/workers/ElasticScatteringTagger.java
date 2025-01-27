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
public class ElasticScatteringTagger extends DataWorker {

  Schema recparticleSchema = null;

  public ElasticScatteringTagger(){
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
  }

  @Override
  public boolean init(DataSource src) {
      //throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
      return true;
  }

  public void cleanArr(double[] arr, double length){
      for(int i=0;i<length;i++){arr[i]=0;}
  }

  public void fillRECPart(Bank RECPart, int pindex, double[] part){
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
  }

  @Override
  public void execute(DataEvent event) {
      Bank bpart = new Bank(recparticleSchema,4096);
      ((Event) event).read(bpart);

      int nVarsPart=9;
      double[] p=new double[nVarsPart];
      double[] el=new double[nVarsPart];

      for(int row=0;row<bpart.getRows();row++){
        for (int row2 = 0; row2 < bpart.getRows(); row2++) {


          cleanArr(p,nVarsPart);
          cleanArr(el,nVarsPart);

          fillRECPart(bpart,row2,p);
          fillRECPart(bpart,row,el);


          if( Math.abs(p[4])>4000 && p[5]==1){
            if(Math.abs(el[4])>=2000 && Math.abs(el[4])<4000 && el[5]==-1 ){

              double phiDif=p[3]-el[3];
              if(phiDif>170 && phiDif<190 ){
                ((Event) event).setEventTag(11);

              }//good phi dif

            }//good el candidate
          }//good recoil nucleon candidate

        }//loop over nucleon candidates
      }//loop over electron candidates
    }

}