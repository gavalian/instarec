/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.validation;

import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Bank;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.instarec.utils.DataEntry;
import j4np.instarec.utils.DataList;
import j4np.utils.io.OptionParser;
import twig.data.BarChartBuilder;
import twig.data.DataGroup;
import twig.data.H1F;
import twig.data.H2F;
import twig.data.TDirectory;
import twig.graphics.TGCanvas;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

/**
 *
 * @author tyson
 */
public class DenoisingValidator {
    
    public DenoisingValidator(){
    }

    //run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar
    //read plots in j4shell with eg TwigStudio.browser("plots/OnePion.twig")

    public void process(String file, int limEvs,String endName, double threshold){
      //output training sample size per sector
      int count = 0;
      double[] all= new double[7];
      double[] all_1=new double[7];
      double[] all_0=new double[7];
      double[] acc= new double[7];
      double[] acc_1=new double[7];
      double[] acc_0=new double[7];

      HipoReader r = new HipoReader();
      Event ev = new Event();
      r.open(file);

      Leaf dchits = new Leaf(27, 1, "i", 4800);
      //don't count if set limEvs=-1
      if(limEvs==-1){count=-2;}

      H1F hall = new H1F("Sector (All)", 6,0.5,6.5);
      hall.attr().setLineColor(1);
      hall.attr().setLineWidth(3);
      hall.attr().setTitleX("Sector");

      H1F hall_0 = new H1F("Sector (All, BG)", 6,0.5,6.5);
      hall_0.attr().setLineColor(1);
      hall_0.attr().setLineWidth(3);
      hall_0.attr().setTitleX("Sector");

      H1F hall_1 = new H1F("Sector (All, Signal)", 6,0.5,6.5);
      hall_1.attr().setLineColor(1);
      hall_1.attr().setLineWidth(3);
      hall_1.attr().setTitleX("Sector");

      H1F hacc = new H1F("Sector (Accepted)", 6,0.5,6.5);
      hacc.attr().setLineColor(2);
      hacc.attr().setLineWidth(3);
      hacc.attr().setTitleX("Sector");

      H1F hacc_0 = new H1F("Sector (Accepted, BG)", 6,0.5,6.5);
      hacc_0.attr().setLineColor(2);
      hacc_0.attr().setLineWidth(3);
      hacc_0.attr().setTitleX("Sector");

      H1F hacc_1 = new H1F("Sector (Accepted, Signal)", 6,0.5,6.5);
      hacc_1.attr().setLineColor(2);
      hacc_1.attr().setLineWidth(3);
      hacc_1.attr().setTitleX("Sector");
    
      while(r.hasNext() && count<limEvs){
        //don't count if set limEvs=-1
        if(limEvs!=-1){
          count++;
        }

        r.nextEvent(ev);
        ev.read(dchits);
        for (int row=0;row<dchits.getRows();row++){
          int detector_type = dchits.getInt(0, row); //no get byte function??
          int sect = dchits.getInt(1, row); 
          int order = dchits.getInt(4, row); 
          double prob = dchits.getDouble(5, row); 
          if(detector_type==6){
           all[0]=all[0]+1;
           all[sect]=all[sect]+1;
           hall.fill(sect);
           if(prob>threshold){
             acc[0]=acc[0]+1;
             acc[sect]=acc[sect]+1;
             hacc.fill(sect);
           }

           if(order==0){
             all_0[0]=all_0[0]+1;
             all_0[sect]=all_0[sect]+1;
             hall_0.fill(sect);
             if(prob>threshold){
               acc_0[0]=acc_0[0]+1;
               acc_0[sect]=acc_0[sect]+1;
               hacc_0.fill(sect);
             }
           } else{
             all_1[0]=all_1[0]+1;
             all_1[sect]=all_1[sect]+1;
             hall_1.fill(sect);
             if(prob>threshold){
               acc_1[0]=acc_1[0]+1;
               acc_1[sect]=acc_1[sect]+1;
               hacc_1.fill(sect);
             }
           }

          }//if rows are from dc
        }//rows in bank
      }//while reading

      double[] ratio= new double[7];
      double[] ratio_0= new double[7];
      double[] ratio_1= new double[7];
      for(int i=0;i<7;i++){
        ratio[i]=acc[i]/all[i];
        ratio_0[i]=acc_0[i]/all_0[i];
        ratio_1[i]=acc_1[i]/all_1[i];
      }

      System.out.printf("\n\n All Hits %.3f BG %.3f Signal %.3f\n",all[0],all_0[0],all_1[0]);
      System.out.printf(" Acc Hits %.3f BG %.3f Signal %.3f\n",acc[0],acc_0[0],acc_1[0]);
      System.out.printf(" Ratio      %.3f BG %.3f Signal %.3f\n",ratio[0],ratio_0[0],ratio_1[0]);

      String strratio=String.format("%.3f",ratio[0]);
      String strratio_0=String.format("%.3f",ratio_0[0]);
      String strratio_1=String.format("%.3f",ratio_1[0]);

      TDirectory.export("plots/DenoiserValidation"+endName+".twig","/ai/",hall);
      TDirectory.export("plots/DenoiserValidation"+endName+".twig","/ai/",hall_0);
      TDirectory.export("plots/DenoiserValidation"+endName+".twig","/ai/",hall_1);
      TDirectory.export("plots/DenoiserValidation"+endName+".twig","/ai/",hacc);
      TDirectory.export("plots/DenoiserValidation"+endName+".twig","/ai/",hacc_0);
      TDirectory.export("plots/DenoiserValidation"+endName+".twig","/ai/",hacc_1);

      BarChartBuilder b = new BarChartBuilder();
      b.addEntry("All Hits e^- ",all[0],all_0[0],all_1[0]);
      b.addEntry("Accepted Hits (Ratio: "+strratio+"%, "+strratio_0+"%, "+strratio_1+"% )",acc[0],acc_0[0],acc_1[0]);
      b.setTitleY("Counts");
      b.setColors(new int[]{1,2});
      b.setLabels(new String[]{"All Hits","BG Hits","Signal Hits"});
      DataGroup b2 = b.build();
      TDirectory.export("plots/DenoiserValidation"+endName+".twig","/ai/",b2);
      TGCanvas c = new TGCanvas(1000,1000);
      //for(DataSet ds : group.getData()) c.draw(ds, "same");
      c.view().region().draw(b2);//.showLegend(0.05, 0.95);
      c.view().region().showLegend(0.05, 0.95);
      c.repaint();
      
    }
    
    public static void main(String[] args){
        
      System.out.println("\n\n----- starting denoising validator ");

      double threshold=0.001;

      String fName="/Users/tyson/data_repo/denoising_data/extracted_back_denoised_0f_03f.h5";
      String endName="_03";

      DenoisingValidator dp = new DenoisingValidator();

      dp.process(fName,-1,endName,threshold);

      fName="/Users/tyson/data_repo/denoising_data/extracted_back_denoised_0f_24f.h5";
      endName="_24";
      dp.process(fName,-1,endName,threshold);

    }
  }