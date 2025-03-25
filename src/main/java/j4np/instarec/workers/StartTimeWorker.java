package j4np.instarec.workers;

import j4np.data.base.DataActor;
import j4np.data.base.DataEvent;
import j4np.data.base.DataSource;
import j4np.data.base.DataWorker;
import j4np.hipo5.data.Bank;
import j4np.hipo5.data.Event;
import j4np.hipo5.data.Leaf;
import j4np.hipo5.data.Schema;
import j4np.hipo5.data.Schema.SchemaBuilder;
import j4np.instarec.utils.NeuralModel;
import twig.data.H1F;
import twig.data.TDirectory;
import java.util.ArrayList;

/**
 * StartTimeWorker class that calculates and saves a histogram of start times.
 */
public class StartTimeWorker extends DataWorker {

    @Override
    public boolean init(DataSource src) {
        return true;
    }

    private ArrayList<Double> startTimes = new ArrayList<>();
    private static final String OUTPUT_TWIG_FILE = "start_time_distribution.twig"; // Output file

    public void CalculateStartTime(Leaf partout, Leaf Startout, short row) {
        double c = 29.9792; // Speed of light in cm/ns
        double t_flight = partout.getDouble(28, row);
        double path = partout.getDouble(27, row);
        double t_0 = t_flight - (path / c);
        System.out.println("Start Time: " + t_0);

        Startout.putDouble(0, 0, t_0);
        startTimes.add(t_0);
    }

    public short CalculateHighestMomentumElectron(Leaf partout, Leaf Startout, Leaf trackbank) {
        double maxP = 0;
        short maxPRow = -1;

        for (short i = 0; i < partout.getRows(); i++) {
            if (partout.getInt(1, i) == 11 && partout.getInt(4, i) == -1) {
                double px = partout.getDouble(6, i);
                double py = partout.getDouble(7, i);
                double pz = partout.getDouble(8, i);
                double rho = (double) Math.sqrt(px * px + py * py + pz * pz);
                short pindex = partout.getShort(0, i);
                double status = trackbank.getInt(0, pindex);
                if (rho > maxP && status>0 ) {
                    maxP = rho;
                    maxPRow = i;
                }
            }
        }

        if (maxPRow != -1) {
            System.out.println("Highest momentum electron found in row: " + maxPRow);
        } else {
            System.out.println("  ");
        }
        return maxPRow;
    }

    @Override
    public void execute(DataEvent event) {
        Leaf partout = new Leaf(32200, 1, "sifssf3f3f6f9f3ffff", 4096);
        ((Event) event).read(partout);

        Leaf Startout = new Leaf(32200, 4, "d", 4096);
        Startout.setRows(1);

        Leaf trackbank = new Leaf(32000,1,"i",4096);
        ((Event) event).read(trackbank);

        short maxPRow = CalculateHighestMomentumElectron(partout, Startout, trackbank);
        if (maxPRow != -1) {
            CalculateStartTime(partout, Startout, maxPRow);
            //System.out.println("StartTimes array size: " + startTimes.size());
        } else {
            Startout.putDouble(0, 0, -999);
        }

        ((Event) event).write(Startout);

       //System.out.println("🔹 DEBUG: finish() method called.");
        //System.out.println("Total startTimes size: " + startTimes.size());

        //if (startTimes.isEmpty()) {
            //System.out.println("No start times recorded, skipping histogram save.");
            //return;
        //}

        // Create histogram
        //H1F hist = new H1F("Start Time Distribution", 1000, 0.0, 10000.0);
        //for (double t : startTimes) {
            //hist.fill(t);
        //}
        //String filePath = System.getProperty("user.dir") + "/" + OUTPUT_TWIG_FILE;
        //System.out.println("🔹 DEBUG: Saving histogram to " + filePath);
        // Save histogram using `TDirectory.export()`
        //TDirectory.export(OUTPUT_TWIG_FILE, "/plots/StartTime", hist);
        //System.out.println("Histogram saved as " + OUTPUT_TWIG_FILE);
    }
}


