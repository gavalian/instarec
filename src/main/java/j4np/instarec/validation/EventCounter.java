/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.validation;

import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.utils.io.OptionParser;
import twig.data.GraphErrors;
import twig.data.H1F;
import twig.data.H2F;
import twig.data.StatNumber;
import twig.data.TDirectory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.IOException;

//for priting
import j4np.hipo5.data.Bank;

/**
 *
 * @author tyson
 */
public class EventCounter {
    
    public EventCounter(){
    }

    public static int[] convertL1Trigger(long bits) {
      int[] trigger = new int[32];
  
      // System.out.printf("%X - %X\n", bits,bits&0xF);
      for (int i = 0; i < trigger.length; i++) {
        trigger[i] = 0;
        if (((bits >> i) & (1L)) != 0L)
          trigger[i] = 1;
        // System.out.println(Arrays.toString(trigger));
      }
      return trigger;
    }

    public int hasTriggerEl(Bank triggerBank){
      long bits = triggerBank.getLong("trigger", 0);
      int[] L1trigger = convertL1Trigger(bits);
      return L1trigger[0];
    }

    public int hasPredEl(Leaf pred_part, double th){
      int hasEl=0;
      for(int row=0;row<pred_part.getRows();row++){
        short charge = pred_part.getShort(4,row);
        float resp=(float)pred_part.getDouble(2,row);
        int pid=pred_part.getInt(1,row);
        //if(resp>th && charge==-1){
        if(pid==11){
          hasEl=1;
        }
      }
      return hasEl;
    }

    public int hasRECEl(Bank RECPart){
      int hasEl=0;
      for(int row=0;row<RECPart.getRows();row++){
        int pid=RECPart.getInt("pid", row);
        //if(resp>th && charge==-1){
        if(pid==11){
          hasEl=1;
        }
      }
      return hasEl;
    }


    public void process(String file, int limEvs,double threshold){

      //output training sample size per sector
      int count = 0;
      float nRECEl=0, nInstaEl=0, nTrig=0;

      HipoReader r = new HipoReader();
      Event ev = new Event();
      r.open(file);
      

      /*for printing */
      Bank recpart = r.getBank("REC::Particle");
      Bank triggerbank = r.getBank("RUN::config");
      Leaf pred_part = new Leaf(32, 3, "i", 1200);

      //don't count if set limEvs=-1
      if(limEvs==-1){count=-2;}
      
      while(r.hasNext() && count<limEvs){
        //don't count if set limEvs=-1
        if(limEvs!=-1){
          count++;
        }

        r.nextEvent(ev);
        //Leaf pred_part = ev.readLeaf(1,12,32,3);
        ev.read(pred_part,32,3);
        ev.read(recpart);
        ev.read(triggerbank);

        nTrig+=hasTriggerEl(triggerbank);
        nRECEl+=hasRECEl(recpart);
        nInstaEl+=hasPredEl(pred_part, threshold);
      }

      float ratio_pred=nInstaEl/nTrig;
      float ratio_rec=nRECEl/nTrig;
      float ratio_predrec=nInstaEl/nRECEl;

      System.out.printf("\n\n Nb Trig %f, Nb REC %f, Nb Insta %f \n",nTrig,nRECEl,nInstaEl);
      System.out.printf("ratio pred %f, ratio rec %f, ratio pred to rec %f \n",ratio_pred,ratio_rec,ratio_predrec);


    }

    public static void main(String[] args){
        
      System.out.println("\n\n----- starting event counter ");
      String fName="wRadPhotons.h5";

      double resp_threshold=0.025; 
        
      EventCounter dp = new EventCounter();

      dp.process(fName,-1,resp_threshold);


     
        
    }

}