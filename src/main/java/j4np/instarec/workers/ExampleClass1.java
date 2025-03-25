import j4np.hipo5.data.Bank;
import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.hipo5.data.Leaf
import java.lang.Math;
import twig.data.TDirectory;
import twig.data.GraphErrors;
import twig.data.H1F;
import twig.graphics.TGCanvas;

public class ExampleClass1 {

  public ExampleClass1(){}

  public static void main(String[] args) {

    // Variable to hold the highest momentum
    int pid;
    float StartTime;
    short maxPRow = -1;


    String filePath2 = "w.h5";

    H1F htmaster = new H1F("Event Start time from Data", 500, 0, 200);
    htmaster.attr().setLineColor(4);
    //ht.attr().setFillColor(2);
    htmaster.attr().setLineWidth(1);
    htmaster.attr().setTitleX("T from data [ns]");

    H2F hPaddleTime = new H2F("hPaddleTime", 100, -40000, 40000, 100, 60, 110);
    //hPaddle.attr().setLineColor(4);
    //ht.attr().setFillColor(2);
    //hPaddle.attr().setLineWidth(1);
    hPaddleTime.attr().setTitleX("Paddle");
    hPaddleTime.attr().setTitleY("Time");

    H2F hMomentumTime = new H2F("hMomentumTime", 500, 0, 10, 500, 60, 110);
    //hPaddle.attr().setLineColor(4);
    //ht.attr().setFillColor(2);
    //hPaddle.attr().setLineWidth(1);
    hMomentumTime.attr().setTitleX("Momentum");
    hMomentumTime.attr().setTitleY("Time");

    H1F ht = new H1F("Event Start time (Momentum below 4 GeV)", 250, 80, 110);
    ht.attr().setLineColor(2);
    //ht.attr().setFillColor(2);
    ht.attr().setLineWidth(1);
    ht.attr().setTitleX("T from data [ns]");

    H1F ht2 = new H1F("Event Start time (Momentum between 4-7 GeV)", 250, 80, 110);
    ht2.attr().setLineColor(3);
    //ht2.attr().setFillColor(2);
    ht2.attr().setLineWidth(1);
    ht2.attr().setTitleX("T from data [ns]");

    H1F ht3 = new H1F("Event Start time (Momentum above 7 GeV)", 250, 80, 110);
    ht3.attr().setLineColor(4);
    //ht3.attr().setFillColor(2);
    ht3.attr().setLineWidth(1);
    ht3.attr().setTitleX("T from data [ns]");

    H1F htpredmaster = new H1F("Predicted Event Start time", 500, 0, 200);
    htpredmaster.attr().setLineColor(2);
    //htpred.attr().setFillColor(2);
    htpredmaster.attr().setLineWidth(1);
    htpredmaster.attr().setTitleX("T pred [ns]");

    H1F htpredmaster2 = new H1F("Predicted Event Start time", 500, 0, 200);
    htpredmaster2.attr().setLineColor(3);
    //htpred.attr().setFillColor(2);
    htpredmaster2.attr().setLineWidth(1);
    htpredmaster2.attr().setTitleX("T pred [ns]");

    H1F htpredmaster3 = new H1F("Predicted Event Start time", 500, 0, 200);
    htpredmaster3.attr().setLineColor(4);
    //htpred.attr().setFillColor(2);
    htpredmaster3.attr().setLineWidth(1);
    htpredmaster3.attr().setTitleX("T pred [ns]");

    H1F htpredmaster4 = new H1F("Predicted Event Start time", 500, 0, 200);
    htpredmaster4.attr().setLineColor(5);
    //htpred.attr().setFillColor(2);
    htpredmaster4.attr().setLineWidth(1);
    htpredmaster4.attr().setTitleX("T pred [ns]");

    H1F htpredmaster5 = new H1F("Predicted Event Start time", 500, 0, 200);
    htpredmaster5.attr().setLineColor(6);
    //htpred.attr().setFillColor(2);
    htpredmaster5.attr().setLineWidth(1);
    htpredmaster5.attr().setTitleX("T pred [ns]");

    H1F htpredmaster6 = new H1F("Predicted Event Start time", 500, 0, 200);
    htpredmaster6.attr().setLineColor(7);
    //htpred.attr().setFillColor(2);
    htpredmaster6.attr().setLineWidth(1);
    htpredmaster6.attr().setTitleX("T pred [ns]");

    H1F ht1pred1 = new H1F("ht1pred", 250, 80, 110);
    ht1pred1.attr().setLineColor(2);
    //htpred.attr().setFillColor(2);
    ht1pred1.attr().setLineWidth(1);
    ht1pred1.attr().setTitleX("T pred [ns]");

    H1F ht1pred2 = new H1F("ht1pred2", 250, 80, 110);
    ht1pred2.attr().setLineColor(8);
    //htpred2.attr().setFillColor(2);
    ht1pred2.attr().setLineWidth(1);
    ht1pred2.attr().setTitleX("T pred [ns]");

    H1F ht1pred3 = new H1F("ht1pred3", 250, 80, 110);
    ht1pred3.attr().setLineColor(4);
    //htpred3.attr().setFillColor(2);
    ht1pred3.attr().setLineWidth(1);
    ht1pred3.attr().setTitleX("T pred [ns]");

    H1F ht2pred1 = new H1F("ht2pred1", 250, 80, 110);
    ht2pred1.attr().setLineColor(2);
    //htpred3.attr().setFillColor(2);
    ht2pred1.attr().setLineWidth(1);
    ht2pred1.attr().setTitleX("T pred [ns]");

    H1F ht2pred2 = new H1F("ht2pred2", 250, 80, 110);
    ht2pred2.attr().setLineColor(8);
    //htpred3.attr().setFillColor(2);
    ht2pred2.attr().setLineWidth(1);
    ht2pred2.attr().setTitleX("T pred [ns]");

    H1F ht2pred3 = new H1F("ht2pred3", 250, 80, 110);
    ht2pred3.attr().setLineColor(4);
    //htpred3.attr().setFillColor(2);
    ht2pred3.attr().setLineWidth(1);
    ht2pred3.attr().setTitleX("T pred [ns]");

    H1F ht3pred1 = new H1F("ht3pred1", 250, 80, 110);
    ht3pred1.attr().setLineColor(2);
    //htpred3.attr().setFillColor(2);
    ht3pred1.attr().setLineWidth(1);
    ht3pred1.attr().setTitleX("T pred [ns]");

    H1F ht3pred2 = new H1F("ht3pred2", 250, 80, 110);
    ht3pred2.attr().setLineColor(8);
    //htpred3.attr().setFillColor(2);
    ht3pred2.attr().setLineWidth(1);
    ht3pred2.attr().setTitleX("T pred [ns]");

    H1F ht3pred3 = new H1F("ht3pred3", 250, 80, 110);
    ht3pred3.attr().setLineColor(4);
    //htpred3.attr().setFillColor(2);
    ht3pred3.attr().setLineWidth(1);
    ht3pred3.attr().setTitleX("T pred [ns]");

    H1F ht4pred1 = new H1F("ht4pred1", 250, 80, 110);
    ht4pred1.attr().setLineColor(2);
    //htpred3.attr().setFillColor(2);
    ht4pred1.attr().setLineWidth(1);
    ht4pred1.attr().setTitleX("T pred [ns]");

    H1F ht4pred2 = new H1F("ht4pred2", 250, 80, 110);
    ht4pred2.attr().setLineColor(8);
    //htpred3.attr().setFillColor(2);
    ht4pred2.attr().setLineWidth(1);
    ht4pred2.attr().setTitleX("T pred [ns]");

    H1F ht4pred3 = new H1F("ht4pred3", 250, 80, 110);
    ht4pred3.attr().setLineColor(4);
    //htpred3.attr().setFillColor(2);
    ht4pred3.attr().setLineWidth(1);
    ht4pred3.attr().setTitleX("T pred [ns]");

    H1F ht5pred1 = new H1F("ht5pred1", 250, 80, 110);
    ht5pred1.attr().setLineColor(2);
    //htpred3.attr().setFillColor(2);
    ht5pred1.attr().setLineWidth(1);
    ht5pred1.attr().setTitleX("T pred [ns]");

    H1F ht5pred2 = new H1F("ht5pred2", 250, 80, 110);
    ht5pred2.attr().setLineColor(8);
    //htpred3.attr().setFillColor(2);
    ht5pred2.attr().setLineWidth(1);
    ht5pred2.attr().setTitleX("T pred [ns]");

    H1F ht5pred3 = new H1F("ht5pred3", 250, 80, 110);
    ht5pred3.attr().setLineColor(4);
    //htpred3.attr().setFillColor(2);
    ht5pred3.attr().setLineWidth(1);
    ht5pred3.attr().setTitleX("T pred [ns]");

    H1F ht6pred1 = new H1F("ht6pred1", 250, 80, 110);
    ht6pred1.attr().setLineColor(2);
    //htpred3.attr().setFillColor(2);
    ht6pred1.attr().setLineWidth(1);
    ht6pred1.attr().setTitleX("T pred [ns]");

    H1F ht6pred2 = new H1F("ht6pred2", 250, 80, 110);
    ht6pred2.attr().setLineColor(8);
    //htpred3.attr().setFillColor(2);
    ht6pred2.attr().setLineWidth(1);
    ht6pred2.attr().setTitleX("T pred [ns]");

    H1F ht6pred3 = new H1F("ht6pred3", 250, 80, 110);
    ht6pred3.attr().setLineColor(4);
    //htpred3.attr().setFillColor(2);
    ht6pred3.attr().setLineWidth(1);
    ht6pred3.attr().setTitleX("T pred [ns]");

    H1F hmomP = new H1F("First Particle in REC::Particle Momentum (REC bank)", 1000, 0, 10);
    hmomP.attr().setLineColor(2);
    hmomP.attr().setFillColor(2);
    hmomP.attr().setLineWidth(3);
    hmomP.attr().setTitleX("P [GeV]");

    H1F hmom = new H1F("Highest Electron Momentum (partout leaf)", 1000, 0, 10);
    hmom.attr().setLineColor(2);
    hmom.attr().setFillColor(2);
    hmom.attr().setLineWidth(3);
    hmom.attr().setTitleX("P [GeV]");

    H1F htdiff = new H1F("Difference between Start times", 100, 0, 200);
    htdiff.attr().setLineColor(2);
    htdiff.attr().setFillColor(2);
    htdiff.attr().setLineWidth(3);
    htdiff.attr().setTitleX("T predicted [ns]");

    HipoReader r = new HipoReader(filePath2);

    Event e = new Event();
    //show the list of schemas
    // r.getSchemaFactory().show();
    Bank[] banks = r.getBanks("REC::Particle", "REC::Event", "REC::Track");
    Leaf partout = new Leaf(32200, 1, "sifssf3f3f6f9f3ffff", 4096);
    Leaf Startout = new Leaf(32200, 4, "d", 4096);
    Leaf trackbank = new Leaf(32000,1,"i",4096);
    Leaf FTOFclusters = new Leaf(32200,3,"i",4096);
    int nEvs=0;
    //double maxP = 0;
    while (r.hasNext()) {
      r.nextEvent(e);
      e.read(banks);
      e.read(Startout, 32200, 4);
      e.read(partout, 32200, 1);
      e.read(trackbank, 32000, 1);
      e.read(FTOFclusters, 32200, 3);

      double maxP = 0;
      int electron_no=0;
      short sector = -1;
      short component = 0;
      short status = 0;
      short pindex = -1;
      //show contents of rec particle
      if(nEvs==0){
        banks[0].show();
        banks[1].show();
      }

      // If there is an electron, fill the histogram with the event start time.

      for (short i = 0; i < banks[0].getRows(); i++) {
        pid = banks[0].getInt("pid", i);
        if (pid == 11) {
          electron_no += 1;
        }
      }

      if(electron_no>0){
        float px = banks[0].getFloat("px", 0);
        float py = banks[0].getFloat("py", 0);
        float pz = banks[0].getFloat("pz", 0);
        float momentum = (float) Math.sqrt(px * px + py * py + pz * pz);

        //short component = banks[2].getInt("component", j);
        //double status =  banks[2].getInt("status", j);
        //double sector = banks[2].getInt("sector", j);
        //System.out.println("Momentum of first particle in event " + momentum);
        hmomP.fill(momentum);
        for (short j = 0; j < banks[1].getRows(); j++) {
            // Get momentum information (using px, py, pz for 3D momentum)
            //System.out.println("Momentum in part bank" + momentum);
            StartTime = banks[1].getFloat("startTime", j);
            htmaster.fill(StartTime);
            if (momentum < 4) {
              ht.fill(StartTime);
            } else if (momentum > 4 && momentum < 7) {
              ht2.fill(StartTime);
            } else {
              ht3.fill(StartTime);
            }

        }

      }

      for (short i = 0; i < partout.getRows(); i++) {
        if (partout.getInt(1, i) == 11 && partout.getInt(4, i) == -1) {
          double px = partout.getDouble(6, i);
          double py = partout.getDouble(7, i);
          double pz = partout.getDouble(8, i);
          //Calculate absolute value of 3-momentum:
          double rho = (double) Math.sqrt(px * px + py * py + pz * pz);

          pindex = partout.getShort(0, i);
          status = trackbank.getShort(0, pindex);
          component = FTOFclusters.getShort(3, pindex);
          sector = FTOFclusters.getShort(1, pindex);
          //Update electron with max momentum dynamically:
          if (rho > maxP && status>0) {
            maxP = rho;
            //System.out.println("Max Electron momentum = " + maxP);
          }
        }
      }
      if (maxP > 0) {
        hmom.fill((float) maxP);  // Add the momentum to the histogram
        //hPaddle.fill((float) component);
      }
      for (short i = 0; i < Startout.getRows(); i++) {
        //System.out.println("maxP = " + maxP);

        double t = Startout.getDouble(0, i);

        if (t>0 && maxP>0 && status>0 && component!=0) {
          hPaddleTime.fill(component, t);

          hMomentumTime.fill(maxP, t);
        }
        if (t>0 && maxP>0 && status>0 && sector==1) {
          htpredmaster.fill(t);
          //System.out.println("Max El momentum: " + maxP);
          if(maxP<4){
            //System.out.println("Start Time: " + t);
            ht1pred1.fill(t);}
          else if(maxP>4 && maxP<7){
            ht1pred2.fill(t);}
          else{ht1pred3.fill(t);}

        }

        if (t>0 && maxP>0 && status>0 && sector==2) {
          htpredmaster2.fill(t);
          if(maxP<4){
            //System.out.println("Start Time: " + t);
            ht2pred1.fill(t);}
          else if(maxP>4 && maxP<7){
            ht2pred2.fill(t);}
          else{ht2pred3.fill(t);}
        }

        if (t>0 && maxP>0 && status>0 && sector==3) {
          htpredmaster3.fill(t);
          if(maxP<4){
            //System.out.println("Start Time: " + t);
            ht3pred1.fill(t);}
          else if(maxP>4 && maxP<7){
            ht3pred2.fill(t);}
          else{ht3pred3.fill(t);}
        }

        if (t>0 && maxP>0 && status>0 && sector==4) {
          htpredmaster4.fill(t);
          if(maxP<4){
            //System.out.println("Start Time: " + t);
            ht4pred1.fill(t);}
          else if(maxP>4 && maxP<7){
            ht4pred2.fill(t);}
          else{ht4pred3.fill(t);}
        }

        if (t>0 && maxP>0 && status>0 && sector==5) {
          htpredmaster5.fill(t);
          if(maxP<4){
            //System.out.println("Start Time: " + t);
            ht5pred1.fill(t);}
          else if(maxP>4 && maxP<7){
            ht5pred2.fill(t);}
          else{ht5pred3.fill(t);}
        }

        if (t>0 && maxP>0 && status>0 && sector==6) {
          htpredmaster6.fill(t);
          if(maxP<4){
            //System.out.println("Start Time: " + t);
            ht6pred1.fill(t);}
          else if(maxP>4 && maxP<7){
            ht6pred2.fill(t);}
          else{ht6pred3.fill(t);}
        }
      }
      nEvs++;
    }
    double entries = hmom.getEntries();
    //double entries2 = htpred.getEntries();
    //double entries3 = htpred2.getEntries();
    //double entries4 = htpred3.getEntries();
    double entries5 = hmomP.getEntries();
    double entries6 = ht.getEntries();
    double entries7 = ht2.getEntries();
    double entries8 = ht3.getEntries();



    /*
    TGCanvas canvas = new TGCanvas();
    canvas.setTitle("Multiple Histograms");
    canvas.view().region().draw(htpredmaster);
    canvas.repaint();


    TGCanvas canvas = new TGCanvas();
    canvas.view().divide(3,2);
    canvas.setTitle("Multiple Histograms");
    canvas.view().cd(0).region().draw(htpredmaster);
    canvas.view().cd(1).region().draw(htpredmaster2);
    canvas.view().cd(2).region().draw(htpredmaster3);
    canvas.view().cd(3).region().draw(htpredmaster4);
    canvas.view().cd(4).region().draw(htpredmaster5);
    canvas.view().cd(5).region().draw(htpredmaster6);
    htpredmaster.attr().setLegend("Sector 1");
    htpredmaster2.attr().setLegend("Sector 2");
    htpredmaster3.attr().setLegend("Sector 3");
    htpredmaster4.attr().setLegend("Sector 4");
    htpredmaster5.attr().setLegend("Sector 5");
    htpredmaster6.attr().setLegend("Sector 6");
    canvas.cd(0).view().region().showLegend(0.5,0.98);
    canvas.cd(1).view().region().showLegend(0.5,0.98);
    canvas.cd(2).view().region().showLegend(0.5,0.98);
    canvas.cd(3).view().region().showLegend(0.5,0.98);
    canvas.cd(4).view().region().showLegend(0.5,0.98);
    canvas.cd(5).view().region().showLegend(0.5,0.98);
    canvas.repaint();

     */

    ht1pred1.attr().setLegend("Sector 1 P<4");
    ht1pred2.attr().setLegend("Sector 1 4<P<7");
    ht1pred3.attr().setLegend("Sector 1 P>7");

    TGCanvas canvas = new TGCanvas();
    canvas.view().divide(3,2);
    canvas.setTitle("Multiple Histograms");
    canvas.cd(0).view().region().draw(ht1pred1).draw(ht1pred2,"same").draw(ht1pred3,"same");
    canvas.cd(0).view().region().showLegend(0.5,0.98);
    canvas.cd(0).repaint();

    ht2pred1.attr().setLegend("Sector 2 P<4");
    ht2pred2.attr().setLegend("Sector 2 4<P<7");
    ht2pred3.attr().setLegend("Sector 2 P>7");

    canvas.cd(1).view().region().draw(ht2pred1).draw(ht2pred2,"same").draw(ht2pred3,"same");
    canvas.cd(1).view().region().showLegend(0.,0.98);
    canvas.cd(1).repaint();

    ht3pred1.attr().setLegend("Sector 3 P<4");
    ht3pred2.attr().setLegend("Sector 3 4<P<7");
    ht3pred3.attr().setLegend("Sector 3 P>7");

    canvas.cd(2).view().region().draw(ht3pred1).draw(ht3pred2,"same").draw(ht3pred3,"same");
    canvas.cd(2).view().region().showLegend(0,0.98);
    canvas.cd(2).repaint();

    ht4pred1.attr().setLegend("Sector 4 P<4");
    ht4pred2.attr().setLegend("Sector 4 4<P<7");
    ht4pred3.attr().setLegend("Sector 4 P>7");

    canvas.cd(3).view().region().draw(ht4pred1).draw(ht4pred2,"same").draw(ht4pred3,"same");
    canvas.cd(3).view().region().showLegend(0.,0.98);
    canvas.cd(3).repaint();

    ht5pred1.attr().setLegend("Sector 5 P<4");
    ht5pred2.attr().setLegend("Sector 5 4<P<7");
    ht5pred3.attr().setLegend("Sector 5 P>7");

    canvas.cd(4).view().region().draw(ht5pred1).draw(ht5pred2,"same").draw(ht5pred3,"same");
    canvas.cd(4).view().region().showLegend(0.,0.98);
    canvas.cd(4).repaint();

    ht6pred1.attr().setLegend("Sector 6 P<4");
    ht6pred2.attr().setLegend("Sector 6 4<P<7");
    ht6pred3.attr().setLegend("Sector 6 P>7");

    canvas.cd(5).view().region().draw(ht6pred1).draw(ht6pred2,"same").draw(ht6pred3,"same");
    canvas.cd(5).view().region().showLegend(0,0.98);
    canvas.cd(5).repaint();
    /*
    ht.attr().setLegend("P<4 (data)");
    ht2.attr().setLegend("4<P<7 (data)");
    ht3.attr().setLegend("P>7 (data)");

    canvas.cd(1).view().region().draw(ht).draw(ht2,"same").draw(ht3,"same");
    canvas.cd(1).view().region().showLegend(0.5,0.98);
    canvas.cd(1).repaint();


    TGCanvas canvas2 = new TGCanvas();
    canvas2.setTitle("Multiple Histograms");
    htmaster.attr().setLegend("Data");
    htpredmaster.attr().setLegend("Prediction");
    canvas2.view().region().draw(htpredmaster,"same").draw(htmaster,"same");
    canvas2.view().region().showLegend(0.5,0.98);
    canvas2.repaint();
    */

    //System.out.println("entries in mom:" + entries);
    //System.out.println("entries in mom from REC Bank:" + entries5);
    //System.out.println("entries in tpred:" + entries2 + entries3 + entries4);
    //System.out.println("entries in t:" + entries6 + entries7 + entries8);
    // Save histogram using `TDirectory.export()`
    TDirectory.export("start_time_Exp1.twig", "/plots/StartTime", ht);
    TDirectory.export("start_time_Exp2.twig", "/plots/StartTime", ht2);
    TDirectory.export("start_time_Exp3.twig", "/plots/StartTime",	ht3);
    //TDirectory.export("start_time_Pred1.twig", "/plots/StartTime", htpred);
    //TDirectory.export("start_time_Pred2.twig", "/plots/StartTime", htpred2);
    //TDirectory.export("start_time_Pred3.twig", "/plots/StartTime", htpred3);
    TDirectory.export("start_time_MomentumExp.twig", "/plots/StartTime", hmomP);
    TDirectory.export("start_time_MomentumPred.twig", "/plots/StartTime", hmom);
    TDirectory.export("start_time_MomentumPred.twig", "/plots/StartTime", hmom);
    System.out.println("Histogram saved as start_time_New.twig");
    System.out.println("Histogram saved as start_time_Exp1.twig");

  }
}
