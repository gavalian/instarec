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

/**
 * StartTimeWorker class that calculates and saves a histogram of start times.
 */
public class StartTimeWorker extends DataWorker {

    @Override
    public boolean init(DataSource src) {
        return true;
    }


    public void CalculateStartTime(Leaf partout, Leaf Startout, short row) {
        // Speed of light in cm/ns
        double c = 29.9792;
        // Get time of flight and path from partout bank written by ElPID worker:
        double t_flight = partout.getFloat(27, row);
        double path = partout.getFloat(28, row);
        // Calculate start time for the event assuming beta=1:
        double t_0 = t_flight - (path / c);
        System.out.println("Start Time: " + t_0);

        Startout.putDouble(0, 0, t_0);
    }

    public short CalculateHighestMomentumElectron(Leaf partout, Leaf Startout) {
        // Finds the highest momentum electron in any given event.
        float maxP = 0;
        short maxPRow = -1;

        for (short i = 0; i < partout.getRows(); i++) {
            if (partout.getInt(1, i) == 11 && partout.getInt(4, i) == -1) {
                float px = partout.getFloat(5, i);
                float py = partout.getFloat(6, i);
                float pz = partout.getFloat(7, i);
                //Calculate absolute value of 3-momentum:
                float rho = (float) Math.sqrt(px * px + py * py + pz * pz);
                //Update electron with max momentum dynamically:
                if (rho > maxP) {
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

        Leaf Startout = new Leaf(32200, 4, "d", 4096);
        Startout.setRows(1);

        short maxPRow = CalculateHighestMomentumElectron(partout, Startout);
        if (maxPRow != -1) {
            // If there is at least one electron in the event:
            CalculateStartTime(partout, Startout, maxPRow);
            System.out.println("StartTimes array size: " + startTimes.size());
        } else {
            // If no electron, put the start time as -999.
            Startout.putDouble(0, 0, -999);
        }
        // Write the event.
        ((Event) event).write(Startout);
    }
    }


