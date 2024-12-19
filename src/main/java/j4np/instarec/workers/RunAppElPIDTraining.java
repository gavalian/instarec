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
import j4np.instarec.data.ElPIDDataProvider;
import j4np.instarec.networks.TrainingElPID;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author tyson
 */
public class RunAppElPIDTraining {
    public static DataFrame createFrames(int count){
        DataFrame<Event>  frame = new DataFrame<>();        
        for(int i = 0; i < count; i++) frame.addEvent(new Event());
        return frame;
    }
    
    public static List<DataActor>  createActors(int nactors, int nframes, List<DataWorker> workers){
        List<DataActor> actors = new ArrayList<>();
        for(int a = 0; a < nactors; a++){
            DataActor actor = new DataActor();
            DataFrame frame = RunAppElPIDTraining.createFrames(nframes);
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

        String pathToClusterFinder = "etc/networks/clusterfinder/cf";
        
        HipoWriter w = HipoWriter.create("w.h5", r);
        
        DataActorStream stream = new DataActorStream();
        
        stream.setSource(r).setSync(w);
        
        ConverterWorker   convert = new ConverterWorker();
        DriftChamberWorker  dcwrk = new DriftChamberWorker();
        TrackFinderWorker  finder = new TrackFinderWorker();
        ClusterFinderWorkerECAL  ecalfinder = new ClusterFinderWorkerECAL(pathToClusterFinder);
        ClusterFinderWorkerFTOF  ftoffinder = new ClusterFinderWorkerFTOF(pathToClusterFinder);
        ConverterWorkerHTCC  htcc = new ConverterWorkerHTCC();
        ConverterWorkerParticleCFTraining   convertParticleCFTraining = new ConverterWorkerParticleCFTraining();
        
        finder.initNetworks();
        
        List<DataWorker>  workers = Arrays.asList(convert,dcwrk,finder,ecalfinder,ftoffinder,htcc,convertParticleCFTraining);
        
        List<DataActor>   actors = RunAppElPIDTraining.createActors(1, 128, workers);
        
        stream.addActor(actors);//.addActor(convert2);//.addActor(convert3).addActor(convert4);
        
        stream.run();

        ElPIDDataProvider dp = new ElPIDDataProvider();
        dp.process("w.h5", "training_data/ElPIDTrain",150000);

        double desiredThreshold=0.99;

        System.out.println("\n\n\nTraining v1");
        String dataPath = "training_data/ElPIDTrain";
        String networkPath = "etc/networks/ElPID/ElPID";
        TrainingElPID.trainNetwork(dataPath, networkPath);
        TrainingElPID.testNetwork(dataPath, networkPath, 0,desiredThreshold);
        for (int j = 1; j < 7; j++) {
          System.out.printf("\nTransfer training for sector %d\n", j);
          TrainingElPID.transferTraining(dataPath, networkPath, j);
          TrainingElPID.testNetwork(dataPath, networkPath, j,desiredThreshold);
        }
        
    }
}
