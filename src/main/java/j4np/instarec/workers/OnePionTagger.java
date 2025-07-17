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
public class OnePionTagger extends DataWorker {

  Schema recparticleSchema = null;
  double beamE=10.6;
  double targetMass=0.938272;

  public OnePionTagger(double be,double m){
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
    beamE=be;
    targetMass=m;
  }

  @Override
  public boolean init(DataSource src) {
      //throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
      return true;
  }

  public double square(double a){
    return a*a;
  }


  public void calcExc(double[] pip, double[] el, double[] exc,double beamE,double tm){
      
    double elE=Math.sqrt(square(el[1])+square(0.000511));
    double pipE=Math.sqrt(square(pip[1])+square(0.13957));
    double pM=tm;

    double IM = Math.sqrt(square(pipE+elE)- ( square(pip[6]+el[6]) + square(pip[7]+el[7]) + square(pip[8]+el[8]) ));
    double px_m = -1.0*(el[6]+pip[6]);
    double py_m = -1.0*(el[7]+pip[7]);
    double pz_m = beamE - (el[8]+pip[8]);
    double p_m=Math.sqrt(square(px_m) + square(py_m)+square(pz_m));
    double pxp_m = px_m/p_m;
    double pyp_m = py_m/p_m;
    double E_m = beamE + tm - (elE+pipE);
    double MM2 = square(E_m) - (square(px_m) + square(py_m)+square(pz_m));

    exc[0]=IM;
    exc[1]=Math.sqrt(MM2);
    exc[2]=Math.sqrt(pxp_m*pxp_m + pyp_m*pyp_m);
    exc[3]=(MM2-pM*pM)/(2*pM);
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

      int nVarsPart=9, nVarsExc=4;
      double[] pip=new double[nVarsPart];
      double[] el=new double[nVarsPart];
      double[] exc=new double[nVarsExc];

      for(int row=0;row<bpart.getRows();row++){
        for (int row2 = 0; row2 < bpart.getRows(); row2++) {


          cleanArr(pip,nVarsPart);
          cleanArr(el,nVarsPart);
          //cleanArr(exc,nVarsExc);

          fillRECPart(bpart,row2,pip);
          fillRECPart(bpart,row,el);


          if(Math.abs(pip[4])>=2000 && Math.abs(pip[4])<4000 && pip[0]==211){
            if(Math.abs(el[4])>=2000 && Math.abs(el[4])<4000 && el[5]==-1 ){
              
              //calcExc(pip, el, exc,beamE,targetMass);
              //cut on missing mass, large for plots and fits
              if(exc[1]>0.2 && exc[1]<2.1){
                ((Event) event).setEventTag(211);
              }

            }//good el candidate
          }//good pi+ candidate

        }//loop over pi+ candidates
      }//loop over electron candidates
    }
}