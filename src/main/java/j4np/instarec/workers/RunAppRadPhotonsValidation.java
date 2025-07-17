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
import j4np.instarec.validation.TwoPionValidator;
import j4np.instarec.validation.TwoPionFTelValidator;
import j4np.instarec.validation.KLambdaValidator;
import j4np.instarec.validation.TwoPionNoElectronValidator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author tyson
 */
public class RunAppRadPhotonsValidation {
    public static DataFrame createFrames(int count){
        DataFrame<Event>  frame = new DataFrame<>();        
        for(int i = 0; i < count; i++) frame.addEvent(new Event());
        return frame;
    }
    
    public static List<DataActor>  createActors(int nactors, int nframes, List<DataWorker> workers){
        List<DataActor> actors = new ArrayList<>();
        for(int a = 0; a < nactors; a++){
            DataActor actor = new DataActor();
            DataFrame frame = RunAppRadPhotonsValidation.createFrames(nframes);
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
        //String file="/Users/tyson/data_repo/trigger_data/rga/run_valid_2.h5";
        String file="/Users/tyson/data_repo/trigger_data/rga/run_validTwoPion_noPimID_inbending.h5";
        
        
        for(int i=1;i<3;i++){
          
        
          String outName="wvalid_TwoPion_noPimID.h5"; //TwoPion KLambda
          if(i==1){
            outName="wvalid_RadPhotons_TwoPion_noPimID.h5";
          } else if (i==2){
            file="/Users/tyson/data_repo/trigger_data/rga/run_validTwoPionNoEl_inbending.h5";
            outName="wvalid_RadPhotons_TwoPionNoEl_noPimID.h5";
          }
          
          HipoReader r = new HipoReader(file);
          HipoWriter w = HipoWriter.create(outName, r); //RadPhotons
          
          DataActorStream stream = new DataActorStream();
          
          stream.setSource(r).setSync(w);

          String pathToClusterFinder = "etc/networks_rga/clusterfinder/cf";
          String pathToElPID = "etc/networks_rga/ElPID/ElPID"; //RadPhotons
          if(i!=0){
            pathToElPID = "etc/networks_rga/ElPID/ElPIDRadPhotons"; //RadPhotons
          }
          double threshold=0.1;        

          ConverterWorker   convert = new ConverterWorker();
          DriftChamberWorker  dcwrk = new DriftChamberWorker();
          TrackFinderWorker  finder = new TrackFinderWorker();
          ClusterFinderWorkerECAL  ecalfinder = new ClusterFinderWorkerECAL(pathToClusterFinder);
          ClusterFinderWorkerFTOF  ftoffinder = new ClusterFinderWorkerFTOF(pathToClusterFinder);
          ConverterWorkerHTCC  htcc = new ConverterWorkerHTCC();
          ElPIDWorker elPID  = new ElPIDWorker(pathToElPID,threshold);
          ConverterWorkerParticleCFTraining   convertParticleCFTraining = new ConverterWorkerParticleCFTraining();
          //TwoPionTagger   TwoPionTagger = new TwoPionTagger(10.6,0.938272);
          //OnePionTagger   OnePionTagger = new OnePionTagger(10.6,0.938272);
          
          finder.initNetworks();
          
          List<DataWorker>  workers = Arrays.asList(convert,dcwrk, finder, ecalfinder, ftoffinder, htcc, elPID, convertParticleCFTraining);
          
          List<DataActor>   actors = RunAppRadPhotonsValidation.createActors(4, 128, workers);
          actors.get(0).setBenchmark(1);
          stream.addActor(actors);//.addActor(convert2);//.addActor(convert3).addActor(convert4);
          
          stream.run();

        }

        String fName="wvalid_TwoPion_noPimID.h5";

        String endName="_noPimID_elPsup2_th0p025_matchTrack";//_elPsup2"; //_phiCut5 eg 175-185, 10 otherwise trackChi2l350_vzl20_6SL

        double resp_threshold=0.025; //0.075
        double beamE=10.6;

        TwoPionValidator dp = new TwoPionValidator();
        //TwoPionFTelValidator dp = new TwoPionFTelValidator();
        //TwoPionNoElectronValidator dp = new TwoPionNoElectronValidator();

        //dp.process(fName,-1,endName,resp_threshold,beamE,0,0);
        //dp.process(fName,-1,endName+"_wFid",resp_threshold,beamE,1,0);
        //dp.process(fName,-1,endName+"_wFidTight",resp_threshold,beamE,2,0);

        fName="wvalid_RadPhotons_TwoPion_noPimID.h5";
        endName="_radPhotonTraining"+endName;

        dp.process(fName,-1,endName,resp_threshold,beamE,0,0);
        dp.process(fName,-1,endName+"_wFid",resp_threshold,beamE,1,0);
        dp.process(fName,-1,endName+"_wFidTight",resp_threshold,beamE,2,0);

        String fNameNoEl="wvalid_TwoPionNoEl_noPimID.h5";

        String endNameNoEl="_noPimID_matchTrack";//_elPsup2"; //_phiCut5 eg 175-185, 10 otherwise trackChi2l350_vzl20_6SL

        TwoPionNoElectronValidator dpNoEl = new TwoPionNoElectronValidator();

        //dp.process(fNameNoEl,-1,endNameNoEl,resp_threshold,beamE,0);

        fNameNoEl="wvalid_RadPhotons_TwoPionNoEl_noPimID.h5";
        endNameNoEl="_radPhotonTraining"+endNameNoEl;

        dpNoEl.process(fNameNoEl,-1,endNameNoEl,resp_threshold,beamE,0);

        
        
    }
}
