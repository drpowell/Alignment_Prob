
package fuzzyLZ;

import java.io.*;

import common.*;

class FuzzyDriver implements Serializable {
    int numIterations;
    char alphabet[];
    int DEBUG;
    String fname;
    int imageFreq;
    int checkpointFreq;
    int statsFreq;

    Params p;
    char[] str;

    public static void main(String args[]) {
        CommandLine cmdLine = new CommandLine();
	cmdLine.addInt("iterations", 1, "Number of iterations");
	cmdLine.addString("alphabet", "atgc", "Alphabet used by the sequence.");
	cmdLine.addInt("debug", 2, "Debug level (higher gives more verbose output)");
	cmdLine.addInt("imageSize",  800, "Maximum Image size in pixels");
	cmdLine.addInt("imageFreq", 0, "Save an image every <n> seconds.  (0 - to disable)");
	cmdLine.addInt("checkFreq", 300, "Save a checkpoint every <n> seconds.  (0 - to disable)");
	cmdLine.addInt("statsFreq", 300, "Display some stats every <n> seconds.  (0 - to disable)");

	// These options are for Matches_Sparse
	cmdLine.addInt("hashSize",   20, "Window size to use for constructing hashtable (0 - for full N^2 algorithm)");
	cmdLine.addInt("computeWin", 10,  "Number of cells to activate on a hashtable hit");
	cmdLine.addInt("cutML", 4,  "When (cell_value - base_cell > cutML) then cell is killed. (in bits)");
	cmdLine.addBoolean("plotActive", false, "true: plot only active cells, false: plot cell values");


	cmdLine.addBoolean("resume", false, "Resume from a checkpoint");

	args = cmdLine.parseLine(args);

	if (args==null || args.length != 1) {
	    System.err.println("Usage: java fuzzyLZ.FuzzyLZ [options] <seqFile|checkpointFile>\n" + cmdLine.usage());
	    System.exit(1);
	}

	boolean resume = cmdLine.getBooleanVal("resume");

	FuzzyDriver me = null;

	if (resume) {
	    try {
		System.out.println("Reloading from checkpoint...");
		File f = new File(args[0]);
		ObjectInputStream is = new ObjectInputStream(new FileInputStream(f));
		me = (FuzzyDriver)is.readObject();
		System.out.println("Successfully loaded checkpoint...");
	    } catch (Exception e) {
		System.err.println("Unable to resume from checkpoint: "+e);
		System.exit(1);
	    }
	    // We can re-set some of the parameters here.
	    // Only reset ones that are defined on this commandline, ie. don't use defaults.
	    if (cmdLine.optionSet("iterations"))  me.numIterations   = cmdLine.getIntVal("iterations");
	    if (cmdLine.optionSet("imageFreq"))   me.imageFreq       = cmdLine.getIntVal("imageFreq");
	    if (cmdLine.optionSet("checkFreq"))   me.checkpointFreq  = cmdLine.getIntVal("checkFreq");
	    if (cmdLine.optionSet("statsFreq"))   me.statsFreq       = cmdLine.getIntVal("statsFreq");
	    if (cmdLine.optionSet("imageSize"))   
		FuzzyLZ.img_width  = FuzzyLZ.img_height = cmdLine.getIntVal("imageSize");
	    if (cmdLine.optionSet("debug")) {
		me.DEBUG = cmdLine.getIntVal("debug");
		FuzzyLZ.DEBUG = me.DEBUG;
	    }
	    
	} else {
	    me = new FuzzyDriver();

	    me.numIterations = cmdLine.getIntVal("iterations");
	    me.alphabet = cmdLine.getStringVal("alphabet").toCharArray();
	    me.DEBUG = cmdLine.getIntVal("debug");
	    FuzzyLZ.DEBUG = me.DEBUG;
	    FuzzyLZ.img_width  = FuzzyLZ.img_height = cmdLine.getIntVal("imageSize");
	    me.imageFreq      = cmdLine.getIntVal("imageFreq");
	    me.checkpointFreq = cmdLine.getIntVal("checkFreq");
	    me.statsFreq      = cmdLine.getIntVal("statsFreq");
	    Matches_Sparse.def_winSize    = cmdLine.getIntVal("hashSize");
	    Matches_Sparse.def_computeWin = cmdLine.getIntVal("computeWin");
	    Matches_Sparse.def_cutML      = cmdLine.getIntVal("cutML");
	    Matches_Sparse.def_plotActive = cmdLine.getBooleanVal("plotActive");

	    me.fname = args[0];
	    DNA seq = DNA.guess_format(me.fname);
	    if (seq == null) {
		System.err.println("Unable to read sequence file");
		System.exit(1);
	    }
	    me.str = seq.sequence();
	    Misc.printf("Sequence length = %d\n",me.str.length);

	    me.p = new Params();
	}
	me.go(resume);
    }

    

    FuzzyDriver() { };

    Seq_Model m1, m2;
    double totd1, totd2;

    int iteration;
    int inner_i;

    void init_iteration() {
	m1 = new Markov0_DNA(4);
	m2 = new FuzzyLZ( p, new Markov0_DNA(4), str, 4);
	//m2 = new FuzzyLZ( p, new Uniform(4), str, 4);

	totd1=0;
	totd2=0;
    }

    void inner_loop() {
	char c = str[inner_i];
	double d1 = m1.update(c, inner_i);
	double d2 = m2.update(c, inner_i);
		
	totd1 += d1;
	totd2 += d2;
	Misc.my_assert(d2 > 0, "Bugger! -ve bits to encode char");
		
	if (DEBUG>=3)
	    Misc.printf("%s %03d: m1=%.2f m2=%.2f tot(m1)=%.2f tot(m2)=%.2f\n", 
			new Misc.VarArgs(c).add(inner_i).add(d1).add(d2).add(totd1).add(totd2));
	
    }
    
    void go (boolean resume) {
	if (resume) {
	    System.out.println("Resuming from iteration="+iteration+" at character="+inner_i+" of "+str.length);
	} else {
	    iteration = 0;
	}

	long last_checkpoint = System.currentTimeMillis();
	long last_image      = System.currentTimeMillis();
	long last_stats      = System.currentTimeMillis();

	for ( ; iteration<numIterations; iteration++) {
	    if (DEBUG>=1)
		Misc.printf("\n\nIteration %d\nParams:\n%s\n", new Object[] { new Integer(iteration), p });

	    if (!resume) {
		init_iteration();
		inner_i = 0;
	    }
	    resume = false;
		
	    for ( ; inner_i<str.length; inner_i++) {
		// Display some stats?
		if (statsFreq>0 && (System.currentTimeMillis() - last_stats)/1000 >= statsFreq) {
		    last_stats = System.currentTimeMillis();
		    ((FuzzyLZ)m2).display_stats();
		}

		// Save a checkpoint?
		if (checkpointFreq>0 &&  (System.currentTimeMillis() - last_checkpoint)/1000 >= checkpointFreq) {
		    try {
			File f = new File(new File(fname+"-checkpoint-"+iteration+"-"+inner_i+".obj").getName());
			if (DEBUG>=1) 
			    System.out.println("Saving checkpoint : " + f.getAbsoluteFile());
			ObjectOutputStream o = new ObjectOutputStream(new FileOutputStream(f));
			o.writeObject(this);
		    } catch (IOException e) {
			System.err.println("Failed to save checkpoint: "+e);
		    }
		    last_checkpoint = System.currentTimeMillis();
		}

		// Save an image?
		if (imageFreq>0 && (System.currentTimeMillis() - last_image)/1000 >= imageFreq) {
		    ((FuzzyLZ)m2).plot.save(Misc.sprintf("out%07d.ppm", new Misc.VarArgs(inner_i)), "Seq: " + fname);
		    last_image = System.currentTimeMillis();
		}

		inner_loop();
	    }

	    if (DEBUG>=0)
		Misc.printf("Total for m1 = %.4f  Total for m2 = %.4f\n", totd1, totd2);
	    
	    if (DEBUG>=1 || statsFreq>0)
		((FuzzyLZ)m2).display_stats();	

	    if (DEBUG>=2)
		((FuzzyLZ)m2).display();

	    p = ((FuzzyLZ)m2).counts_to_params();

	    ((FuzzyLZ)m2).plot.save(Misc.sprintf("out-iter%02d.final.ppm", new Misc.VarArgs(iteration)), 
				    "Seq: " + fname);
	}
    }
}
