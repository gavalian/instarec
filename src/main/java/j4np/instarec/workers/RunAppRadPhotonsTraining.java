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
import j4np.instarec.data.RadPhotonDataProvider;
import j4np.instarec.networks.TrainingElPIDRadPhotons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author tyson
 */
public class RunAppRadPhotonsTraining {
    public static DataFrame createFrames(int count){
        DataFrame<Event>  frame = new DataFrame<>();        
        for(int i = 0; i < count; i++) frame.addEvent(new Event());
        return frame;
    }
    
    public static List<DataActor>  createActors(int nactors, int nframes, List<DataWorker> workers){
        List<DataActor> actors = new ArrayList<>();
        for(int a = 0; a < nactors; a++){
            DataActor actor = new DataActor();
            DataFrame frame = RunAppRadPhotonsTraining.createFrames(nframes);
            actor.setWorkes(workers);
            actor.setDataFrame(frame);
            actors.add(actor);
        }
        return actors;
    }
    
    public static void main(String[] args){
        
        //String file = "/Users/gavalian/Work/DataSpace/decoded/clas_006595.evio.00625-00629_DC.hipo";
        //  String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_3_wAIBanks.h5";
        // String file = "/Users/tyson/data_repo/trigger_data/rgd/018777/sorted_out_skim_018777_noPID.hipo";
        //String file = "/Users/tyson/data_repo/trigger_data/sims/claspyth_train/clasdis_62.hipo";
        // String file = "/w/work/clas12/tyson/data_repo/caos/rga/run_test_5407.h5";
        // String file = "/work/clas12/jnp/instarec/irec_005197.evio.00011.h5";
        //String file = "/work/clas12/jnp/instarec/irec_rec_005197.evio.h5";
       // String file="/Users/tyson/data_repo/trigger_data/rga/irec_rec_005197.evio.h5";
        //String file="/Users/tyson/data_repo/trigger_data/rga/run_test_5407.h5";
        
        String outName="wvalid.h5";
        String outNameNeg="wtrainNeg.h5";
        
        for(int i=2;i<2;i++){

          String file="/Users/tyson/data_repo/trigger_data/rga/run_valid_inbending.h5";
          if(i==1){file="/Users/tyson/data_repo/trigger_data/rga/run_trainTwoPion_noPimID_5197.hipo";}
          HipoReader r = new HipoReader(file);
          
          String ot=outName;
          if(i==1){ot=outNameNeg;}
          HipoWriter w = HipoWriter.create(ot, r); 

          DataActorStream stream = new DataActorStream();

          stream.setSource(r).setSync(w);

          String pathToClusterFinder = "etc/networks_rga/clusterfinder/cf";
          String pathToElPID = "etc/networks_rga/ElPID/ElPID";
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

          List<DataActor>   actors = RunAppRadPhotonsTraining.createActors(4, 128, workers);
          actors.get(0).setBenchmark(1);
          stream.addActor(actors);//.addActor(convert2);//.addActor(convert3).addActor(convert4);

          stream.run();
        }

        //data preparation
        String dataPath = "training_data/ElPIDTrainRadPhotons";
        String dataPath_noRad = "training_data/ElPIDTrain"; 

        RadPhotonDataProvider dp = new RadPhotonDataProvider();

        //more negatives than positives
        //int nNegatives=dp.negativesFromTwoPion(outNameNeg,100000,"",dataPath,10.6,true);
        //int nPositives=dp.process(outName,nNegatives,0.05,"",dataPath,false);
        //int nPositives=100001;
        //dp.copyOverNegatives(nPositives,dataPath_noRad,dataPath);

        //some rad photons, some usual positives and negatives
        int nPositives=dp.process(outName,50000,0.05,"",dataPath,true);
        int nPositives2=dp.copyOverPositives(2*nPositives,dataPath_noRad,dataPath);
        int nNegatives=dp.copyOverNegatives(nPositives+nPositives2,dataPath_noRad,dataPath);

        //training
        double desiredThreshold=0.99;

        System.out.println("\n\n\nTraining RadPhotons");
       
        String networkPath = "etc/networks_rga/ElPID/ElPID";
        String networkPathOut = "etc/networks_rga/ElPID/ElPIDRadPhotons";

        //transfer training
        /*for (int j = 1; j < 7; j++) {
          System.out.printf("\nTransfer training for sector %d\n", j);
          TrainingElPIDRadPhotons.transferTraining(dataPath, networkPath,networkPathOut, j);
          TrainingElPIDRadPhotons.testNetwork(dataPath, networkPathOut, j,desiredThreshold);
        }*/

        //train from scratch
        TrainingElPIDRadPhotons.trainNetwork(dataPath,networkPathOut);
        TrainingElPIDRadPhotons.testNetwork(dataPath, networkPathOut, 0,desiredThreshold);
        
        
    }
}
