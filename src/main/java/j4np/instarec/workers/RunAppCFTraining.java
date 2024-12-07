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
import j4np.instarec.data.ClusterMatchingDataProvider;
import j4np.instarec.networks.TrainingClusterFinder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author gavalian
 */
public class RunAppCFTraining {
    public static DataFrame createFrames(int count){
        DataFrame<Event>  frame = new DataFrame<>();        
        for(int i = 0; i < count; i++) frame.addEvent(new Event());
        return frame;
    }
    
    public static List<DataActor>  createActors(int nactors, int nframes, List<DataWorker> workers){
        List<DataActor> actors = new ArrayList<>();
        for(int a = 0; a < nactors; a++){
            DataActor actor = new DataActor();
            DataFrame frame = RunAppCFTraining.createFrames(nframes);
            actor.setWorkes(workers);
            actor.setDataFrame(frame);
            actors.add(actor);
        }
        return actors;
    }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar
    
    public static void main(String[] args){
        
        //String file = "/Users/gavalian/Work/DataSpace/decoded/clas_006595.evio.00625-00629_DC.hipo";
        String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_1_wAIBanks.h5";
        //String file = "/Users/tyson/data_repo/trigger_data/sims/claspyth_train/clasdis_62.hipo";
        HipoReader r = new HipoReader(file);
        
        HipoWriter w = HipoWriter.create("w.h5", r);
        
        DataActorStream stream = new DataActorStream();
        
        stream.setSource(r).setSync(w);
        
        ConverterWorker   convert = new ConverterWorker();
        DriftChamberWorker  dcwrk = new DriftChamberWorker();
        TrackFinderWorker  finder = new TrackFinderWorker();

        ConverterWorkerParticleCFTraining   convertParticleCFTraining = new ConverterWorkerParticleCFTraining();
        
        finder.initNetworks();
        
        List<DataWorker>  workers = Arrays.asList(convert,dcwrk, finder,convertParticleCFTraining);
        
        List<DataActor>   actors = RunAppCFTraining.createActors(1, 128, workers);
        
        stream.addActor(actors);//.addActor(convert2);//.addActor(convert3).addActor(convert4);
        
        stream.run();

        ClusterMatchingDataProvider dp = new ClusterMatchingDataProvider();
        dp.process("w.h5", "training_data/cfTrain",150000);

        String[] charge = new String[2];
        charge[0]="negatives";
        charge[1]="positives";

        for(int i=0;i<2;i++){
          System.out.println("\n\n\nTraining for "+charge[i]);
          String dataPath="training_data/cfTrain_"+charge[i];
          String networkPath="etc/networks/clusterfinder/cf_"+charge[i];
          TrainingClusterFinder.trainNetwork(dataPath,networkPath);
          TrainingClusterFinder.testNetwork(dataPath,networkPath,0,charge[i]);
          for(int j=1;j<7;j++){
            System.out.printf("\nTransfer training for sector %d\n",j);
            TrainingClusterFinder.transferTraining(dataPath,networkPath,j);
            TrainingClusterFinder.testNetwork(dataPath,networkPath,j,charge[i]);
          }
        }
        
    }
}
