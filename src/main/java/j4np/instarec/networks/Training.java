/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.networks;

import deepnetts.net.FeedForwardNetwork;
import j4np.instarec.utils.DataEntry;
import j4np.instarec.utils.DataList;
import j4np.instarec.utils.DeepNettsIO;
import j4np.instarec.utils.DeepNettsTrainer;
import j4np.instarec.utils.EntryTransformer;
import j4np.instarec.utils.NeuralModel;
import java.util.Arrays;
import java.util.Random;
import javax.visrec.ml.data.DataSet;


/**
 *
 * @author gavalian
 */
public class Training {
    
    
    public static DataList generate(int count, int nIn, int nOut){
        DataList data = new DataList();
        Random r = new Random();
        for(int i = 0; i < count; i++){
            float[]  input = new float[nIn];
            float[] output = new float[nOut];
            for(int j = 0 ; j < nIn; j++) input[j] = r.nextFloat()*112.0f;
            int which = r.nextInt(nOut); output[which] = 1.0f;
            data.add(new DataEntry(input,output));
        }
        return data;
    }
    
    public static void trainNetwork(){
        
        FeedForwardNetwork network = DeepNettsTrainer.createClassifier(new int[]{6,24,12,6,3});
        EntryTransformer     trans = new EntryTransformer();
        
        trans.input().add(1, 112); // Add individually 
        trans.input().add(1, 112);
        trans.input().add(1, 112);
        trans.input().add(1, 112);
        trans.input().add(1, 112);
        trans.input().add(1, 112);
        
        trans.output().add(3, 0, 1); // add three of the same
        
        DataList data = Training.generate(25000, 6, 3);
        data.show();
        
        DeepNettsTrainer   trainer = new DeepNettsTrainer(network);
        trainer.updateLayers();        
        DataSet dset = trainer.convert(data, trans);        
        trainer.train(dset, 32);
        trainer.save("classifier.json", trans);
    }
    
    public static void testNetwork(){
        NeuralModel model = NeuralModel.jsonFile("classifier.json");
        System.out.println(model.summary());
        
        DataList data = Training.generate(15, 6, 3);
        
        float[] predicted = new float[3];
        
        for(int j = 0; j < data.getList().size(); j++){
            model.predict(data.getList().get(j).features(), predicted);
            data.getList().get(j).setInfered(predicted);
            
            data.getList().get(j).show();
        }
    }
    
    public static void transferTraining(){
        FeedForwardNetwork network = DeepNettsIO.read("classifier.json");//DeepNettsTrainer.createClassifier(new int[]{12,24,12,6,3});        
        EntryTransformer   transformer = DeepNettsIO.getTransformer("classifier.json");
        
        DeepNettsTrainer   trainer = new DeepNettsTrainer(network);
        trainer.updateLayers();
        
        DataList data = Training.generate(25000, 6, 3);
        data.show();
        
        DataSet dset = trainer.convert(data, transformer);        
        trainer.train(dset, 32);
        trainer.save("classifier_transfer.json", transformer);
    }
    
    public static void main(String[] args){
        
        Training.trainNetwork();
        Training.testNetwork();
        Training.transferTraining();
    }
}
