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

public class RunAppHadronPID {
    public static DataFrame createFrames(int count) {
        DataFrame<Event> frame = new DataFrame<>();
        for (int i = 0; i < count; i++) frame.addEvent(new Event());
        return frame;
    }

    public static List<DataActor> createActors(int nactors, int nframes, List<DataWorker> workers) {
        List<DataActor> actors = new ArrayList<>();
        for (int a = 0; a < nactors; a++) {
            DataActor actor = new DataActor();
            DataFrame frame = RunApp.createFrames(nframes);
            actor.setWorkes(workers);
            actor.setDataFrame(frame);
            actors.add(actor);
        }
        return actors;
    }

    public static void main(String[] args) {
        String file = "/work/clas12/tyson/data_repo/caos/rgd/018326/run_18326_3.h5";
        HipoReader r = new HipoReader(file);

        HipoWriter w = HipoWriter.create("w.h5", r);

        DataActorStream stream = new DataActorStream();

        stream.setSource(r).setSync(w);

        String pathToClusterFinder = "etc/networks/clusterfinder/cf";
        String pathToElPID = "etc/networks/ElPID/ElPID";
        double threshold = 0.1;

        ConverterWorker convert = new ConverterWorker();
        DriftChamberWorker dcwrk = new DriftChamberWorker();
        TrackFinderWorker finder = new TrackFinderWorker();
        ClusterFinderWorkerECAL ecalfinder = new ClusterFinderWorkerECAL(pathToClusterFinder);
        ClusterFinderWorkerFTOF ftoffinder = new ClusterFinderWorkerFTOF(pathToClusterFinder);
        ElPIDWorker elPID = new ElPIDWorker(pathToElPID, threshold);
        StartTimeWorker startTimeWorker = new StartTimeWorker(); 

        finder.initNetworks();

       
        List<DataWorker> workers = Arrays.asList(
                convert,
                dcwrk,
                finder,
                ecalfinder,
                ftoffinder,
                elPID,
                startTimeWorker 
        );

        List<DataActor> actors = RunApp.createActors(1, 128, workers);

        stream.addActor(actors);

        stream.run();
    }
}
