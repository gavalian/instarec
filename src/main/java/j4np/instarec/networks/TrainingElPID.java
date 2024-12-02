/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.networks;

import deepnetts.data.TabularDataSet;
import deepnetts.net.FeedForwardNetwork;
import j4np.instarec.utils.DataEntry;
import j4np.instarec.utils.DataList;
import j4np.instarec.utils.DeepNettsIO;
import j4np.instarec.utils.DeepNettsTrainer;
import j4np.instarec.utils.EJMLModel;
import j4np.instarec.utils.EntryTransformer;
import j4np.instarec.utils.NeuralModel;
import j4np.utils.io.OptionParser;
import j4np.utils.io.TextFileReader;
import twig.data.GraphErrors;
import twig.data.H1F;
import twig.data.StatNumber;
import twig.data.TDirectory;
import twig.graphics.TGCanvas;
import twig.studio.TwigStudio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
import javax.visrec.ml.data.DataSet;

/**
 *
 * @author tyson
 */
public class TrainingElPID {

  static int nPredictedOutput = 2;

  // Will get rid of this, there seems to be a bug in DataList.csv() method
  public static DataList readCSV(String file, int nIn, int nOut) {
    DataList list = new DataList();

    try {
      BufferedReader reader = new BufferedReader(new FileReader(file));
      String line;
      while ((line = reader.readLine()) != null) {

        String[] tokens = line.split(",");
        if (tokens.length != nIn + nOut) {
          throw new IllegalArgumentException("Invalid line length: " + line);
        }

        // Create arrays for nIn and nOut
        float[] inputArray = new float[nIn];
        float[] outputArray = new float[nOut];

        // Fill inputArray
        for (int i = 0; i < nIn; i++) {
          inputArray[i] = Float.parseFloat(tokens[i].trim());
        }

        // Fill outputArray
        for (int i = 0; i < nOut; i++) {
          outputArray[i] = Float.parseFloat(tokens[nIn + i].trim());
        }

        list.add(new DataEntry(inputArray, outputArray));

      }

      reader.close();

    } catch (IOException e) {
      System.err.println("Error reading file: " + e.getMessage());
    } catch (NumberFormatException e) {
      System.err.println("Error parsing number: " + e.getMessage());
    }

    return list;
  }

  public static void trainNetwork(String dataPath, String networkPath) {

    FeedForwardNetwork network = DeepNettsTrainer.createClassifier(new int[] { 35, 35, 20, 10, 5, nPredictedOutput });
    EntryTransformer trans = new EntryTransformer();

    trans.input().add(3, 0, 150000); // calo ADCs
    trans.input().add(9, 1, 500); // ECAL locs
    trans.input().add(9, 0, 10); // ECALl DUs
    trans.input().add(6, 1, 112); // wires
    trans.input().add(8, 0, 35000); // ECALl DUs
    trans.output().add(2, 0, 1); // Electron or not

    // DataList data = DataList.csv(dataPath + ".csv",
    // DataList.range(0, 35), DataList.range(35, 35 + nPredictedOutput));
    // //_outbending

    DataList data = readCSV(dataPath + ".csv", 35, nPredictedOutput);

    data.shuffle();
    data.scan();
    // data.show();

    DeepNettsTrainer trainer = new DeepNettsTrainer(network);
    trainer.updateLayers();
    DataSet dset = trainer.convert(data, trans);
    trainer.train(dset, 1000);
    trainer.save(networkPath + ".json", trans);
  }

  public static void testNetwork(String dataPath, String networkPath, int sector, double effLow) {
    String networkPath_ws = networkPath;
    String dataPath_ws = dataPath;
    if (sector != 0) {
      networkPath_ws = networkPath + "_sector" + String.valueOf(sector);
      dataPath_ws = dataPath + "_sector" + String.valueOf(sector);
    }
    NeuralModel model = NeuralModel.jsonFile(networkPath_ws + ".json");
    // System.out.println(model.summary());

    // DataList data = DataList.csv(dataPath_ws + ".csv",
    // DataList.range(0, 35), DataList.range(35, 35 + nPredictedOutput));
    // //_outbending
    DataList data = readCSV(dataPath_ws + ".csv", 35, nPredictedOutput);
    data.shuffle();

    float[] predicted = new float[nPredictedOutput];

    H1F hRespPos = new H1F("e-", 100, 0, 1);
    hRespPos.attr().setLineColor(2);
    hRespPos.attr().setFillColor(2);
    hRespPos.attr().setLineWidth(3);
    hRespPos.attr().setTitleX("Response");
    H1F hRespNeg = new H1F("Not e-", 100, 0, 1);
    hRespNeg.attr().setLineColor(5);
    hRespNeg.attr().setLineWidth(3);
    hRespNeg.attr().setTitleX("Response");

    for (int j = 0; j < data.getList().size(); j++) {
      model.predict(data.getList().get(j).features(), predicted);
      float[] desired = data.getList().get(j).floatSecond();
      if (j < 10) {
        if (j < 1) {
          System.out.println("Testing -- first 10 entries");
        }
        System.out.println("PID true, pred");
        System.out.printf("%f, %f\n", desired[1], predicted[1]);
      }
      if (desired[1] == 1) {
        hRespPos.fill(predicted[1]);
      } else {
        hRespNeg.fill(predicted[1]);
      }
    }
    TDirectory.export("plots/ElPID" + String.valueOf(sector) + ".twig", "/ai/training/Resp", hRespPos);
    TDirectory.export("plots/ElPID" + String.valueOf(sector) + ".twig", "/ai/training/Resp", hRespNeg);
    TrainingElPID.findBestThreshold(data,model,1,effLow,sector);

  }

  // Labels col 0 is 1 if there's an e-, 0 otherwise
  public static StatNumber[] getMetrics(DataList dle, NeuralModel model, double thresh, int elClass) {

    int nEvents = dle.getList().size();

    int nEls = 0;
    double TP = 0, FP = 0, FN = 0;
    float[] output = new float[nPredictedOutput];
    for (int i = 0; i < nEvents; i++) {
      float[] input = dle.getList().get(i).floatFirst();
      float[] desired = dle.getList().get(i).floatSecond();
      model.predict(dle.getList().get(i).features(), output);
      if (desired[elClass] == 1) {
        nEls++;
        if (output[elClass] > thresh) {
          TP++;
        } else {
          FN++;
        }
      } else {
        if (output[elClass] > thresh) {
          FP++;
        }
      } // Check true label
    }

    StatNumber Pur = new StatNumber(TP, Math.sqrt(TP));
    StatNumber Eff = new StatNumber(TP, Math.sqrt(TP));
    StatNumber FPs = new StatNumber(FP, Math.sqrt(FP));
    StatNumber FNs = new StatNumber(FN, Math.sqrt(FN));
    StatNumber denomPur = new StatNumber(TP, Math.sqrt(TP));
    StatNumber denomEff = new StatNumber(TP, Math.sqrt(TP));
    denomPur.add(FPs);
    denomEff.add(FNs);
    Pur.divide(denomPur);
    Eff.divide(denomEff);
    StatNumber[] mets = new StatNumber[2];
    mets[0] = Pur;
    mets[1] = Eff;
    return mets;
  }

  public static double findBestThreshold(DataList dle, NeuralModel model, int elClass, double effLow, int sector) {

    GraphErrors gEff = new GraphErrors();
    gEff.attr().setMarkerColor(2);
    gEff.attr().setMarkerSize(10);
    gEff.attr().setTitle("Efficiency");
    gEff.attr().setTitleX("Response");
    gEff.attr().setTitleY("Metrics");
    GraphErrors gPur = new GraphErrors();
    gPur.attr().setMarkerColor(5);
    gPur.attr().setMarkerSize(10);
    gPur.attr().setTitle("Purity");
    gPur.attr().setTitleX("Response");
    gPur.attr().setTitleY("Metrics");
    double bestRespTh = 0;
    double bestPuratEffLow = 0;

    // Loop over threshold on the response
    for (double RespTh = 0.01; RespTh < 0.99; RespTh += 0.01) {
      StatNumber metrics[] = getMetrics(dle, model, RespTh, elClass);
      double Pur = metrics[0].number();
      double Eff = metrics[1].number();
      double PurErr = metrics[0].error();
      double EffErr = metrics[1].error();
      gPur.addPoint(RespTh, Pur, 0, PurErr);
      gEff.addPoint(RespTh, Eff, 0, EffErr);
      if (Eff > effLow) {
        if (Pur > bestPuratEffLow) {
          bestPuratEffLow = Pur;
          bestRespTh = RespTh;
        }
      }
    } // Increment threshold on response

    System.out.format("%n Best Purity at Efficiency above %f: %.3f at a threshold on the response of %.3f %n%n",
        effLow, bestPuratEffLow, bestRespTh);

    TDirectory.export("plots/ElPID" + String.valueOf(sector) + ".twig", "/ai/training/Purity", gPur);
    TDirectory.export("plots/ElPID" + String.valueOf(sector) + ".twig", "/ai/training/Efficiency", gEff);

    return bestRespTh;
  }

  public static void transferTraining(String dataPath, String networkPath, int sector) {
    String networkPath_ws = networkPath;
    String dataPath_ws = dataPath;
    if (sector != 0) {
      networkPath_ws = networkPath + "_sector" + String.valueOf(sector);
      dataPath_ws = dataPath + "_sector" + String.valueOf(sector);
    }
    FeedForwardNetwork network = DeepNettsIO.read(networkPath + ".json");
    EntryTransformer transformer = DeepNettsIO.getTransformer(networkPath + ".json");

    DeepNettsTrainer trainer = new DeepNettsTrainer(network);
    trainer.updateLayers();

    // DataList data = DataList.csv(dataPath_ws + ".csv",
    // DataList.range(0, 35), DataList.range(35, 35 + nPredictedOutput));
    // //_outbending
    DataList data = readCSV(dataPath_ws + ".csv", 35, nPredictedOutput);
    data.shuffle();
    data.scan();
    // data.show();

    DataSet dset = trainer.convert(data, transformer);
    trainer.train(dset, 250);
    trainer.save(networkPath_ws + ".json", transformer);
  }

  // run with java -jar target/instarec-1.1.1-jar-with-dependencies.jar -in training_data/ElPIDTrain -name etc/networks/ElPID/ElPID
  // read plots in j4shell with eg TwigStudio.browser("plots/ElPID0.twig");

  public static void main(String[] args) {

    System.out.println("----- starting cluster finder training");
    OptionParser p = new OptionParser();
    p.addRequired("-in", "input name");
    p.addRequired("-name", "network name");
    p.parse(args);

    double desiredThreshold=0.99;

    System.out.println("\n\n\nTraining v1");
    String dataPath = p.getOption("-in").stringValue();
    String networkPath = p.getOption("-name").stringValue();
    TrainingElPID.trainNetwork(dataPath, networkPath);
    TrainingElPID.testNetwork(dataPath, networkPath, 0,desiredThreshold);
    for (int j = 1; j < 7; j++) {
      System.out.printf("\nTransfer training for sector %d\n", j);
      TrainingElPID.transferTraining(dataPath, networkPath, j);
      TrainingElPID.testNetwork(dataPath, networkPath, j,desiredThreshold);
    }

  }
}
