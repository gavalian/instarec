/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.workers;

import j4np.data.base.DataActor;
import j4np.data.base.DataActorStream;
import j4np.data.base.DataFrame;
import j4np.data.base.DataWorker;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.hipo5.io.HipoWriter;
import j4np.instarec.core.DriftChamberWorker;
import j4np.instarec.core.TrackFinderWorker;
import j4np.instarec.validation.ClusterMatchingValidator;
import j4np.instarec.validation.ElPIDValidator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author gavalian
 */
public class RunAppValidation {
    public static DataFrame createFrames(int count){
        DataFrame<Event>  frame = new DataFrame<>();        
        for(int i = 0; i < count; i++) frame.addEvent(new Event());
        return frame;
    }
    
    public static List<DataActor>  createActors(int nactors, int nframes, List<DataWorker> workers){
        List<DataActor> actors = new ArrayList<>();
        for(int a = 0; a < nactors; a++){
            DataActor actor = new DataActor();
            DataFrame frame = RunAppValidation.createFrames(nframes);
            actor.setWorkes(workers);
            actor.setDataFrame(frame);
            actors.add(actor);
        }
        return actors;
    }
    
    public static void main(String[] args){
        
        //String file = "/Users/gavalian/Work/DataSpace/decoded/clas_006595.evio.00625-00629_DC.hipo";
         String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_3_wAIBanks.h5";
        //String file = "/Users/tyson/data_repo/trigger_data/sims/claspyth_train/clasdis_62.hipo";
        HipoReader r = new HipoReader(file);
        
        HipoWriter w = HipoWriter.create("w.h5", r);
        
        DataActorStream stream = new DataActorStream();
        
        stream.setSource(r).setSync(w);

        String pathToClusterFinder = "etc/networks/clusterfinder/cf";
        String pathToElPID = "etc/networks/ElPID/ElPID";
        double threshold=0.1;        
                
        ConverterWorker   convert = new ConverterWorker();
        DriftChamberWorker  dcwrk = new DriftChamberWorker();
        TrackFinderWorker  finder = new TrackFinderWorker();
        ClusterFinderWorkerECAL  ecalfinder = new ClusterFinderWorkerECAL(pathToClusterFinder);
        ClusterFinderWorkerFTOF  ftoffinder = new ClusterFinderWorkerFTOF(pathToClusterFinder);
        ConverterWorkerHTCC  htcc = new ConverterWorkerHTCC();
        ElPIDWorker elPID  = new ElPIDWorker(pathToElPID,threshold);
        ConverterWorkerParticleCFTraining   convertParticleCFTraining = new ConverterWorkerParticleCFTraining();
        
        finder.initNetworks();
        
        List<DataWorker>  workers = Arrays.asList(convert,dcwrk, finder, ecalfinder, ftoffinder, htcc, elPID, convertParticleCFTraining);
        
        List<DataActor>   actors = RunAppValidation.createActors(1, 128, workers);
        actors.get(0).setBenchmark(1);
        stream.addActor(actors);//.addActor(convert2);//.addActor(convert3).addActor(convert4);
        
        stream.run();

        ClusterMatchingValidator valid = new ClusterMatchingValidator();

        String[] chargeSt = new String[2];
        chargeSt[0]="negatives";
        chargeSt[1]="positives";

        //repeating reading of file 3*2*7 times
        // ==> innefficient
        //doing so to avoid writing code for 3*3*2*7 + 2*2*7 histos ...
        for(int lay=0;lay<3;lay++){
          for(int charge=0;charge<2;charge++){
            for(short sector=0;sector<1;sector++){ //7
              valid.process("w.h5",150000,sector,chargeSt[charge],lay);
            }
          }
        }

        float[] respbins = new float[3];
        respbins[0]=(float)0.01;
        respbins[1]=(float)0.01;
        respbins[2]=(float)0.99;

        float[] pbins = new float[3];
        pbins[0]=(float)1;
        pbins[1]=(float)1;
        pbins[2]=(float)9;

        float[] thetabins = new float[3];
        thetabins[0]=(float)5;
        thetabins[1]=(float)5;
        thetabins[2]=(float)35;

        float[] phibins = new float[3];
        phibins[0]=(float)10;
        phibins[1]=(float)-180;
        phibins[2]=(float)180;

        ElPIDValidator validpid = new ElPIDValidator();
        validpid.process("w.h5",150000,0.075,respbins,pbins,thetabins,phibins,(float)0.99,false,"");
        
    }
}
