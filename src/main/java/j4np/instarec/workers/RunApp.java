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
import j4np.utils.io.OptionParser;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

//import j4np.clas12.decoder.Clas12DecoderService;
//import j4np.clas12.decoder.Clas12FitterService;
//import j4np.clas12.decoder.Clas12TranslateService;
import j4np.instarec.core.Hipo2Hipo4Worker;
/**
 *
 * @author gavalian
 */
public class RunApp {
    public static DataFrame createFrames(int count){
        DataFrame<Event>  frame = new DataFrame<>();        
        for(int i = 0; i < count; i++) frame.addEvent(new Event());
        return frame;
    }
    
    public static List<DataActor>  createActors(int nactors, int nframes, List<DataWorker> workers){
        List<DataActor> actors = new ArrayList<>();
        for(int a = 0; a < nactors; a++){
            DataActor actor = new DataActor();
            DataFrame frame = RunApp.createFrames(nframes);
            actor.setWorkes(workers);
            actor.setDataFrame(frame);
            actors.add(actor);
        }
        return actors;
    }
    
    public static void main(String[] args){
        
        OptionParser p = new OptionParser();
        
        p.addOption("-t", "4", "number of threads");
        p.addRequired("-o", "output file");
        
        p.parse(args);
        
        //String file = "/Users/gavalian/Work/DataSpace/decoded/clas_006595.evio.00625-00629_DC.hipo";
        String file = "/Users/tyson/data_repo/trigger_data/rgd/018326/run_18326_3_wAIBanks.h5";
        //String file = "/Users/tyson/data_repo/trigger_data/sims/claspyth_train/clasdis_62.hipo";
        
        String output = p.getOption("-o").stringValue();
        int nt = p.getOption("-t").intValue();
        file = p.getInputList().get(0);
        HipoReader r = new HipoReader(file);
        
        HipoWriter w = HipoWriter.create(output, r);
        
        DataActorStream stream = new DataActorStream();
        
        stream.setSource(r).setSync(w);

        String pathToClusterFinder = "etc/networks_rga/clusterfinder/cf";
        String pathToElPID = "etc/networks_rga/ElPID/ElPID";
        double threshold=0.1;        
        stream.setSource(r).setSync(w);
        
                
        ConverterWorker   convert = new ConverterWorker();
        
        //Clas12DecoderService decoder = new Clas12DecoderService();
        //Clas12FitterService   fitter = new Clas12FitterService();
        //Clas12TranslateService trans = new Clas12TranslateService();
        //decoder.keepEvio = true;
        //trans.setKeepEvio(true);

        DriftChamberWorker  dcwrk = new DriftChamberWorker();
        TrackFinderWorker  finder = new TrackFinderWorker();
        //finder.finder.DO_5_SUPERLAYER = false;
        
        ClusterFinderWorkerECAL  ecalfinder = new ClusterFinderWorkerECAL(pathToClusterFinder);
        ClusterFinderWorkerFTOF  ftoffinder = new ClusterFinderWorkerFTOF(pathToClusterFinder);
        ElPIDWorker elPID  = new ElPIDWorker(pathToElPID,threshold);
        
        finder.initNetworks();
        Hipo2Hipo4Worker     h2h4 = new Hipo2Hipo4Worker(64*1024);
        //List<DataWorker>  workers = Arrays.asList(decoder,fitter,trans,dcwrk, finder,ecalfinder,ftoffinder,elPID,h2h4);
        //List<DataWorker>  workers = Arrays.asList(decoder,fitter,trans,dcwrk, finder,ecalfinder,ftoffinder,elPID,h2h4);
        //List<DataWorker>  workers = Arrays.asList(convert,dcwrk, finder);
        List<DataWorker>  workers = Arrays.asList(convert,dcwrk, finder, ecalfinder, ftoffinder, elPID);
        
        List<DataActor>   actors = RunApp.createActors(nt, 128, workers);
        actors.get(0).setBenchmark(1);
        stream.addActor(actors);//.addActor(convert2);//.addActor(convert3).addActor(convert4);
        
        stream.run();
        
    }
}
