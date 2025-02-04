# instarec
Track Reconstruction with AI

# Usage

The reconstruction works on EVIO files, The code in class RunApp.java processes the file (raw EVIO) and creates a hybrid file with EVIO and H5 containers.

First the EVIO file must be converted to EVIO-6/HIPO format, using

```
prompt> j4np.sh c12 -convert -o input.hipo data.evio
```

Then the reconstruction code must be run on the produced HIPO-5 file by downloading coatjava branch "iss-442":

```
prompt> git clone -b iss-442 https://github.com/JeffersonLab/coatjava.git
prompt> cd coatjava
prompt> ./build-coatjava.sh
prompt> ./coatjava/bin/decoder -i input.hipo -o output.hipo
prompt> ./coatjava/bin/recon-util -i output.hipo -o rec_output.hipo -y etc/services/data-cv.yaml -n 1000
```

At this point the file will contain the hybrid HIPO-5 and reconstrcuted HIPO data.

## Reading the output

Here is an example of how to read InstaRec Banks from the hybrid reconstructed HIPO-5 data:

```
HipoReader r = new HipoReader("infile.h5");
Bank[] b = r.getBanks("TimeBasedTrkg::TBTracks");

Event e = new Event();

while(r.next(e)){
    e.read(b);
    b[0].show();
    Leaf leaf = e.readLeaf(1,12,32000,1);
    leaf.print();
}
```

The readLeaf() method provides the container ids (1,12) which is the same for all H5 leafs, then the leaf ids to read, in the example it reads tracks leaf containing identified tracks by AI.

## Profiling

A very lightweight approach to look at the number of objects per class is Class Histogram.

Just generate class histogram by jmap -histo <PID>.

Use -XX:+PrintClassHistogramBeforeFullGC, -XX:+PrintClassHistogramAfterFullGC to see class histograms when your stop-the-world pause occurs. You can compare the snapshot before/after to see which class instances get collected during the stop-the-world pause.

