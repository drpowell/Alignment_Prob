/*  
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 * 
 * This file is part of FuzzyLZ 
 *
 * FuzzyLZ is a program orginally intended for the
 * compression of DNA sequeces.  It can be viewed as a
 * compression model like Lempel-Ziv 77, but instead of
 * exact matches, allowing matches that contain
 * inserts/deletes/mismatches.
 *
 */


package fuzzyLZ;

import java.util.*;
import java.io.*;

import common.*;

class FuzzyDriver implements Serializable {
    // If any object variables are added here, be use to think about serialisation of them...
    // if the variable is needed after a 'resume', add the variable to the 
    // writeObject() and readObject() functions.

    int numIterations;
    char alphabet[];
    int DEBUG;
    String fname;
    int imageFreq;
    int checkpointFreq;
    int statsFreq;
    String msgFname;

    String fprefix;
    FileWriter msgFile;

    Params p;
    char[] str;
    char[] preStr;
    char[] joinedStr;		// Simple concatenation of preStr + str

    public static void main(String args[]) {
        CommandLine cmdLine = new CommandLine();
	cmdLine.addInt("iterations", 1, "Number of iterations");
	cmdLine.addString("preFile", "", "Prepend 'preFile' to sequence.");
	cmdLine.addString("alphabet", "atgc", "Alphabet used by the sequence.");

	cmdLine.addString("fwdMach", "1state", "Comma separated list of machines to use for forward matches.\n" +
			  "(Use an empty string '' for no forward machines.)\n" +
			  "Supported: 1state,3state");

	cmdLine.addString("revMach", "1state", "Comma separated list of machines to use for reverse matches.\n" +
			  "(Use an empty string '' for no reverse machines.)\n" +
			  "Supported: 1state,3state");

	cmdLine.addInt("debug", 2, "Debug level (higher gives more verbose output)");
	cmdLine.addInt("imageSize",  1024, "Maximum Image size in pixels");
	cmdLine.addInt("imageFreq", 0, "Save an image every <n> seconds.  (0 - to disable)");
	cmdLine.addInt("checkFreq", 300, "Save a checkpoint every <n> seconds.  (0 - to disable)");
	cmdLine.addInt("statsFreq", 300, "Display some stats every <n> seconds.  (0 - to disable)");

	cmdLine.addString("msgFile", "", "Output file for encode length of each character.\n"+
			  "(The default is based on the input file name)");

	// These options are for Matches_Sparse
	cmdLine.addInt("hashSize",   20, "Window size to use for constructing hashtable (0 - for full N^2 algorithm)");
	cmdLine.addInt("computeWin", 10,  "Number of cells to activate on a hashtable hit");
	cmdLine.addInt("cutML", 4,  "When (cell_value - base_cell > cutML) then cell is killed. (in bits)");
	//cmdLine.addBoolean("plotActive", false, "true: plot only active cells, false: plot cell values");


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
	    me.msgFname       = cmdLine.getStringVal("msgFile");
	    Matches_Sparse.def_winSize    = cmdLine.getIntVal("hashSize");
	    Matches_Sparse.def_computeWin = cmdLine.getIntVal("computeWin");
	    Matches_Sparse.def_cutML      = cmdLine.getIntVal("cutML");
	    //Matches_Sparse.def_plotActive = cmdLine.getBooleanVal("plotActive");

	    Vector machs = new Vector();
	    FuzzyLZ.def_numFwd = parseMachineNames(cmdLine.getStringVal("fwdMach"), machs);
	    FuzzyLZ.def_numRev = parseMachineNames(cmdLine.getStringVal("revMach"), machs);
	    FuzzyLZ.MutationModels = (String[])machs.toArray(new String[machs.size()]);

	    String preFile = cmdLine.getStringVal("preFile");
	    if (preFile == "") {
		me.preStr = new char[0];
	    } else {
		DNA seq = DNA.guess_format(preFile);
		if (seq == null) {
		    System.err.println("Unable to read prefix file: '"+preFile+"'");
		    System.exit(1);
		}
		me.preStr = seq.sequence();
		Misc.printf("Pre-sequence length = %d\n",me.preStr.length);
	    }	    
	    
	    me.fname = args[0];
	    DNA seq = DNA.guess_format(me.fname);
	    if (seq == null) {
		System.err.println("Unable to read sequence file");
		System.exit(1);
	    }
	    me.str = seq.sequence();
	    Misc.printf("Sequence length = %d\n",me.str.length);


	    // Display the machines we are going to use
	    for (int i=0; i<FuzzyLZ.def_numFwd; i++)
		System.err.println("FwdMachine["+i+"]: "+FuzzyLZ.MutationModels[i]);
	    for (int i=0; i<FuzzyLZ.def_numRev; i++)
		System.err.println("RevMachine["+i+"]: "+FuzzyLZ.MutationModels[i+FuzzyLZ.def_numFwd]);

	    // Create joinedStr
	    me.joinedStr = new char[me.preStr.length + me.str.length];
	    System.arraycopy(me.preStr, 0, me.joinedStr, 0, me.preStr.length);
	    System.arraycopy(me.str, 0,    me.joinedStr, me.preStr.length, me.str.length);

	    // Get the file prefix to use for output filenames
	    me.fprefix = (new File(me.fname)).getName();
	    me.p = new Params();

	    // Open the file for msglen output
	    try {
		File f;
		if (me.msgFname.equals(""))
		    me.msgFname = me.fprefix + "-msglen.txt";
		f = new File(me.msgFname);
		if (f.exists()) {
		    System.err.println("Output file '"+me.msgFname+"' already exists.");
		    System.exit(1);
		}
		if (!f.createNewFile()) {
		    System.err.println("Unable to create Output file '"+me.msgFname);
		    System.exit(1);
		}
		me.msgFile = new FileWriter(f);
	    } catch (IOException e) {
		System.err.println("Error creating msglen output file: "+e);
		System.exit(1);
	    }
	}
	me.go(resume);
    }

    

    FuzzyDriver() { }

    
    static int parseMachineNames(String l, Vector res) {
	String s[] = l.split(",");
	int num=0;
	for (int i=0; i<s.length; i++) {
	    if (s[i].equals("")) continue;
	    if (s[i].compareToIgnoreCase("1state")==0) {
		res.add("common.Mutation_1State$All");
	    } else if (s[i].compareToIgnoreCase("3state")==0) {
		res.add("common.Mutation_3State$All");
	    } else {
		System.err.println("Unknown machine: '"+s[i]+"'");
		System.exit(1);
	    }
	    num++;
	}
	return num;
    }


    // Note: these object variables must also saved/loaded in the writeObject/readObject functions.

    FuzzyLZ mdl;
    double tot_msglen;

    int iteration;
    int inner_i;

    void init_iteration() {
	// Ensure Params 'p' has no funky numbers
	int numParams = p.get_num();
	for (int i=0; i<numParams; i++) {
	    String name = p.get_name_by_id(i);
	    double v = p.get(name);
	    if (Double.isInfinite(v) || Double.isNaN(v)) {
		System.err.println("Potential parameter problem: '"+name+"' has bad value="+v);
	    }
	}


	mdl = new FuzzyLZ( p, new Markov0_DNA(4), joinedStr, 4, preStr.length);
	//mdl = new FuzzyLZ( p, new Uniform(4), str, 4);

	tot_msglen=0;
    }

    void inner_loop() {
	char c = str[inner_i];
	double d = mdl.update(c, preStr.length + inner_i);

	tot_msglen += d;
	Misc.my_assert(d > 0, "Bugger! -ve bits to encode char");
		
	if (DEBUG>=3)
	    Misc.printf("%s %03d: m=%.2f tot(m)=%.2f\n", 
			new Misc.VarArgs(c).add(inner_i).add(d).add(tot_msglen));

	try {
	    msgFile.write(Misc.sprintf("%s %03d %f\n", new Misc.VarArgs(c).add(inner_i).add(d)));
	    msgFile.flush();
	} catch (IOException e) {
	    System.err.println("Error writing to msglen output file: "+e);
	    System.exit(1);
	}
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
		    System.out.println("Stats as at " + new java.util.Date());
		    mdl.display_stats();
		}

		// Save a checkpoint?
		if (checkpointFreq>0 &&  (System.currentTimeMillis() - last_checkpoint)/1000 >= checkpointFreq) {
		    try {
			String f = fprefix + "-checkpoint-"+iteration+"-"+inner_i+".obj";
			if (DEBUG>=1) 
			    System.out.println("Saving checkpoint : " + f);
			ObjectOutputStream o = new ObjectOutputStream(new FileOutputStream(new File(f)));
			o.writeObject(this);
		    } catch (IOException e) {
			System.err.println("Failed to save checkpoint: "+e);
		    }
		    last_checkpoint = System.currentTimeMillis();
		}

		// Save an image?
		if (imageFreq>0 && (System.currentTimeMillis() - last_image)/1000 >= imageFreq) {
		    mdl.plot.save(Misc.sprintf(fprefix + "-tmp-%02d-%07d.ppm", 
					       new Misc.VarArgs(iteration).add(inner_i)), 
				  "Seq: " + fname);
		    mdl.plotActive.save(Misc.sprintf(fprefix + "-tmpActive-%02d-%07d.ppm", 
						     new Misc.VarArgs(iteration).add(inner_i)), 
					"Seq: " + fname);
		    mdl.plotHits.save(Misc.sprintf(fprefix + "-tmpHits-%02d-%07d.ppm", 
						   new Misc.VarArgs(iteration).add(inner_i)), 
				      "Seq: " + fname);
		    last_image = System.currentTimeMillis();
		}

		inner_loop();
	    }

	    if (DEBUG>=0)
		Misc.printf("Total for mdl = %.4f\n", tot_msglen);
	    
	    if (DEBUG>=1 || statsFreq>0) {
		System.out.println("Stats as at " + new java.util.Date());
		mdl.display_stats();	
	    }

	    if (DEBUG>=2)
		mdl.display();

	    p = mdl.counts_to_params();

	    mdl.plot.save(Misc.sprintf(fprefix + "-final-iter%02d.ppm", 
				       new Misc.VarArgs(iteration)), 
			  "Seq: " + fname);
	    mdl.plotActive.save(Misc.sprintf(fprefix + "-finalActive-iter%02d.ppm", 
					     new Misc.VarArgs(iteration)), 
				"Seq: " + fname);
	    mdl.plotHits.save(Misc.sprintf(fprefix + "-finalHits-iter%02d.ppm", 
					   new Misc.VarArgs(iteration)), 
			      "Seq: " + fname);
	}
    }




    // Write our own serization handler.  
    // Just store everything.  Must re-open the msglen output file
    private void writeObject(java.io.ObjectOutputStream out) throws IOException {
	out.writeInt(numIterations);
	out.writeObject(alphabet);
	out.writeInt(DEBUG);
	out.writeObject(fname);
	out.writeInt(imageFreq);
	out.writeInt(checkpointFreq);
	out.writeInt(statsFreq);
	out.writeObject(msgFname);

	out.writeObject(fprefix);
	out.writeObject(p);
	out.writeObject(str);
	out.writeObject(preStr);


	out.writeObject(mdl);
	out.writeDouble(tot_msglen);
	out.writeInt(iteration);
	out.writeInt(inner_i);
    }

    private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException {
	numIterations  = in.readInt();
	alphabet       = (char[])in.readObject();
	DEBUG          = in.readInt();
	fname          = (String)in.readObject();
	imageFreq      = in.readInt();
	checkpointFreq = in.readInt();
	statsFreq      = in.readInt();
	msgFname       = (String)in.readObject();

	fprefix        = (String)in.readObject();
	p              = (Params)in.readObject();
	str            = (char[])in.readObject();
	preStr         = (char[])in.readObject();

	mdl            = (FuzzyLZ)in.readObject();
	tot_msglen     = in.readDouble();
	iteration      = in.readInt();
	inner_i        = in.readInt();

	// re-Create joinedStr
	joinedStr = new char[preStr.length + str.length];
	System.arraycopy(preStr, 0, joinedStr, 0, preStr.length);
	System.arraycopy(str, 0,    joinedStr, preStr.length, str.length);


	// Open the file for msglen output
	try {
	    File f;
	    f = new File(msgFname);
	    if (!f.exists()) {
		System.err.println("Output file '"+msgFname+"' does not exist. Unable to resume");
		System.exit(1);
	    }
	    msgFile = new FileWriter(f, true);
	    msgFile.write("# Resuming from checkpoint...\n");
	} catch (IOException e) {
	    System.err.println("Error creating msglen output file: "+e);
	    System.exit(1);
	}
    }
}
