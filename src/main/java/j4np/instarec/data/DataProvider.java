/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package j4np.instarec.data;

import j4np.hipo5.data.Event;
import j4np.hipo5.io.HipoReader;
import j4np.hipo5.io.HipoWriter;
import j4np.utils.io.OptionParser;

/**
 *
 * @author gavalian
 */
public class DataProvider {
    
    public DataProvider(){
        
    }
    
    public void process(String file, String output){
        HipoReader r = new HipoReader(file);
        HipoWriter w = HipoWriter.create(output, r);
        Event e = new Event();
        
        while(r.next(e)){
            w.add(e);
        }
        
        w.close();
        
    }
    
    public static void main(String[] args){
        System.out.println("----- starting data provider ");
        OptionParser p = new OptionParser();
        p.addRequired("-o", "output file name");
        p.parse(args);
        
        DataProvider dp = new DataProvider();
        
        dp.process(p.getInputList().get(0), p.getOption("-o").stringValue());
        
    }
}
