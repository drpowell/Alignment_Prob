
package alignCompress;

import java.io.*;
import common.*;

/**
   A model of 2 sequences (for alignments) with counts.

   Characters from each sequence are encoded with a sequence
   specific model.  Matches/changes average the probabilites from these models.

   Uses two parameters: match_cost, change_cost
**/
class Model_SeqAB implements Two_Seq_Model_Counts {
    double match_cost,change_cost;
    Seq_Model modelA, modelB;

    int countIndex;
    final private int matchIndex=0, changeIndex=1;

    public Model_SeqAB(Params p, Seq_Model modelA, Seq_Model modelB, int countIndex) {
        this.countIndex = countIndex;
	this.modelA = modelA;
	this.modelB = modelB;

        if (!p.exists("match_cost")) {
            set_default_costs();
        } else {
            match_cost = p.get("match_cost");
            change_cost = p.get("change_cost");
        }

        normalize_costs();
    }

    public String toString() {
	return this.getClass() + ": match_cost="+match_cost+" change_cost="+change_cost;
    }

    void set_default_costs() {
        match_cost  = -MyMath.log2(9.0);
        change_cost = -MyMath.log2(1.0);
    }

    void normalize_costs() {
        double sum = MyMath.exp2(-match_cost)+MyMath.exp2(-change_cost);
        match_cost  = match_cost + MyMath.log2( sum );
        change_cost = change_cost + MyMath.log2( sum );
    }

    public double encA(char a, int i) {
        return modelA.encodeLen(a, i);
    }

    public double encB(char a, int i) {
        return modelB.encodeLen(a, i);
    }

    public double encBoth(char a, char b, int i, int j) {
	double A_cost = encA(a,i);
	double B_cost = encB(b,j);
	if (a==b) {
	    // Match
	    // Do: P(match) * ( P(char a) + P(char b) ) / 2
	    //System.err.println("enc match = " + ( match_cost + MyMath.logplus(A_cost, B_cost) + 1));
	    return match_cost + MyMath.logplus(A_cost, B_cost) + 1;
	} else {
	    // Change
	    // Do: P(change) * P(char a) * P(char b) * 0.5 * (1/(1-P(char b)) + 1/(1-P(char a)))
	    double aN = MyMath.exp2(-encA(b,i));
	    double bN = MyMath.exp2(-encB(a,j));
	    double norm = -MyMath.log2( 1/(1-aN) + 1/(1-bN) );
	    //System.err.println("enc change = " + (change_cost + A_cost + B_cost + 1 + norm));
	    return change_cost + A_cost + B_cost + 1 + norm;
	}
    }

    public static int required_counts() { return 2; }
    public Params counts_to_params(Counts counts) {
        double sum = counts.counts[countIndex+matchIndex] + 
            counts.counts[countIndex+changeIndex];
        Params par = new Params();
        par.put("match_cost", -MyMath.log2(counts.counts[countIndex+matchIndex]/sum));
        par.put("change_cost", -MyMath.log2(counts.counts[countIndex+changeIndex]/sum));
        return par;
    }

    public void update_count_encA(Counts c, double w, char a, int i) {};
    public void update_count_encB(Counts c, double w, char a, int i) {};
    public void update_count_encBoth(Counts c, double w, char a, char b, int i, int j) {
        if (a==b) {
            c.inc(countIndex+matchIndex, w);
        } else {
            c.inc(countIndex+changeIndex, w);
        }
    }

    public double encode_params(double N) { 
	return Multinomial.MMLparameter_cost(new double[] { MyMath.exp2(-match_cost), 
							    MyMath.exp2(-change_cost) },
					     N);
    }
}

/** BufferModel takes a model and a sequence and the sequence alphabet.
    It precomputes the encoding length of every alphabet character for every 
    position in the sequence.
    The cumulative encoding length is also computed along the sequence.

    This is class probably not practical for large sequences with
    large alphabets.
**/
class BufferModel implements Seq_Model {
    String str;
    char[] alphabet;
    double[][] enc;
    double[] encCumulative;

    BufferModel(Seq_Model model, String str, char[] alpha) {
	this.str = str;
	this.alphabet = alpha;

	enc = new double[str.length()][alphabet.length];
	encCumulative = new double[str.length()+1];
	encCumulative[0] = 0;

	for (int i=0; i<str.length(); i++) {
	    for (int j=0; j<alphabet.length; j++) {
		enc[i][j] = model.encodeLen( alphabet[j], i );
	    }
	    model.update( str.charAt(i), i);
	    encCumulative[i+1] = encCumulative[i] + enc[i][ char2int(str.charAt(i)) ];
	}

	// Test 'model' is efficient and doesn't cheat.
	for (int i=0; i<str.length(); i++) {
	    double s = Double.POSITIVE_INFINITY;
	    for (int j=0; j<alphabet.length; j++) {
		s = MyMath.logplus(s, enc[i][j]);
	    }
	    Misc.assert(Math.abs(s) < 1E-10, "logplus() != 0.  s="+s);
	}
    }

    public String toString() {
	StringBuffer s = new StringBuffer();
	for (int i=0; i<str.length(); i++) {
	    s.append(str.charAt(i) + "    ");
	    for (int j=0; j<alphabet.length; j++) {
		double n = Math.rint(enc[i][j] * 100)/100;
		s.append( alphabet[j] + ": "+n+" ");
	    }
	    s.append("\n");
	}
	return s.toString();
    }

    int char2int(char a) {
	for (int i=0; i<alphabet.length; i++)
	    if (a == alphabet[i])
		return i;
	Misc.assert(false, "Charater '"+a+"' not in defined alphabet");
	return -1;
    }

    public double encodeLen(char a, int i) {
	return enc[i][char2int(a)];
    }

    public double update(char a, int i) {
	return enc[i][char2int(a)];
    }

    public double encodeCumulative(int i) {
	return encCumulative[i];
    }
}

class AlignCompress {
    char[] alphabet;
    String seqA;
    String seqB;

    String paramString;
    int markovOrder;
    int maxIterations;
    int verbose;
    boolean linearCosts;
    boolean localAlign;
    boolean sumAlignments;

    boolean doTraceBack;

    //  static final double PcharEncode = 1.0/170;

    // Encode length l: 0..infinity
    static private double encode_length(double l) {
	Misc.assert(l>=0, "Bad length to encode:"+l);
	return MyMath.logstar_continuous(l+1);
	//return -(l*MyMath.log2(PcharEncode) + MyMath.log2(1-PcharEncode)); // Geometric distribution
    }

    public AlignCompress() {
    }

    // Return the appropriate cell of D[][].  Use i%2 if only keeping 2 rows
    private Mutation_FSM cell(Mutation_FSM D[][], int i, int j) {
	if (doTraceBack) {
	    return D[i][j];
	} else {
	    return D[i%2][j];
	}
    }

    private String getAlignment(Mutation_FSM D[][], 
				Mutation_FSM.TraceBack_Info fcell) {
	Mutation_FSM.TraceBack_Data d = fcell.get_tbdata();
	int i=d.i;
	int j=d.j;

	StringBuffer resA = new StringBuffer();
	StringBuffer resB = new StringBuffer();

	if (i<0 && j<0) {
	    d = fcell.get_from(d);
	    i = d.i;
	    j = d.j;
	}

	String endA = seqA.substring(i);
	String endB = seqB.substring(j);
	
	while (true) {
	    //	    System.err.println("\n\n(i,j)="+i+","+j+"\nA="+resA+"\nB="+resB+"\n");

	    d = ((Mutation_FSM.TraceBack_Info)D[i][j]).get_from(d);

	    if (d==null) break;

	    resA.append((d.i==i) ? '-' : seqA.charAt(d.i));
	    resB.append((d.j==j) ? '-' : seqB.charAt(d.j));

	    i=d.i;
	    j=d.j;
	}

	resA.reverse();
	resB.reverse();

	StringBuffer res = new StringBuffer();

	// Pretty up the alignment.  For local alignments, put the front and end bits on.
	for (int x=0; x< j-(i<j ? i : j); x++) res.append(" ");
	res.append(seqA.substring(0,i));
	if (i>0 || j>0) res.append("   ");
	res.append(resA.toString().toUpperCase());
	res.append("   ");
	res.append(endA);

	res.append("\n");

	for (int x=0; x< i-(i<j ? i : j); x++) res.append(" ");
	res.append(seqB.substring(0,j));
	if (i>0 || j>0) res.append("   ");
	res.append(resB.toString().toUpperCase());
	res.append("   ");
	res.append(endB);

	return res.toString();
    }

    public double doAlign() {
	System.out.println("# SeqA = "+seqA+"\n# SeqB = "+seqB+
			   "\n# markov order="+markovOrder+
			   "\n# max iterations="+maxIterations+
			   "\n# linearCosts="+linearCosts+
			   "\n# sumAlignments="+sumAlignments+
			   "\n# local Alignment="+localAlign+
			   "\n# verbosity="+verbose+
			   "\n# params="+paramString+
			   "\n");

	Params p = new Params();
	p.fromString(paramString);

	BufferModel modelA = new BufferModel( new MarkovN_fitted(markovOrder, alphabet, seqA), 
					      seqA, alphabet );
	BufferModel modelB = new BufferModel( new MarkovN_fitted(markovOrder, alphabet, seqB), 
					      seqB, alphabet );

	//System.err.println("modelA\n"+modelA+"\n");
	//System.err.println("modelB\n"+modelB+"\n");

	int myCounts  = (localAlign ? 1 : 0);
	int mdlCounts = Model_SeqAB.required_counts();

	int fsmCounts = -1;
	if (!linearCosts && !sumAlignments) fsmCounts = Mutation_1State.One.required_counts();
	if (!linearCosts &&  sumAlignments) fsmCounts = Mutation_1State.All.required_counts();
	if ( linearCosts && !sumAlignments) fsmCounts = Mutation_3State.One.required_counts();
	if ( linearCosts &&  sumAlignments) fsmCounts = Mutation_3State.All.required_counts();

	Misc.assert(fsmCounts>=0, "Bad number of fsmCounts");

	int totCounts = myCounts + mdlCounts + fsmCounts;

	double bestDiff = Double.NEGATIVE_INFINITY;
	double lastAlignment = 0;

	int iter=0;
	while (  maxIterations<0 || iter<maxIterations ) {
	    //double PalignChar = 0.9;
	    //if (p.exists("PalignChar")) PalignChar = p.get("PalignChar");

	    int countPos = myCounts;
	    Two_Seq_Model_Counts model = new Model_SeqAB(p, modelA, modelB, countPos);
	    countPos += mdlCounts;
	    Mutation_FSM fsmType=null;

	    if (!linearCosts && !sumAlignments) 
		fsmType = new Mutation_1State.One(model, p, totCounts, countPos);
	    if (!linearCosts &&  sumAlignments) 
		fsmType = new Mutation_1State.All(model, p, totCounts, countPos);
	    if ( linearCosts && !sumAlignments) 
		fsmType = new Mutation_3State.One(model, p, totCounts, countPos);
	    if ( linearCosts &&  sumAlignments) 
		fsmType = new Mutation_3State.All(model, p, totCounts, countPos);

	    Misc.assert(fsmType!=null, "Unable to construct fsmType");


	    countPos += fsmCounts;
	    Misc.assert(countPos == totCounts, "Internal error: countPos!=totCounts");

	    Counts initialCounts = new Counts(totCounts);
	    for (int i=0; i<totCounts; i++) 
		initialCounts.inc(i, 0.5); // Initialise all counts.

	    // Determine if fsmType supports traceBack. ie. Does it implement TraceBack_Info
	    doTraceBack = false;
	    if (fsmType instanceof Mutation_FSM.TraceBack_Info)
		doTraceBack = true;
	    
	    
	    // Setup the DPA matrix
	    Mutation_FSM D[][];
	    Mutation_FSM final_cell;
	    if (!doTraceBack) {
		// No traceback info, only keep 2 rows in D[][]
		D = new Mutation_FSM[2][seqB.length()+1];
		for (int i=0; i<2; i++) {
		    for (int j=0; j<seqB.length()+1; j++) {
			D[i][j] = (Mutation_FSM)fsmType.clone();
		    }
		}
		final_cell = (Mutation_FSM)fsmType.clone();
	    } else {
		// Keep all traceback info, so keep all rows in D[][]
		D = new Mutation_FSM[seqA.length()+1][seqB.length()+1];
		for (int i=0; i<seqA.length()+1; i++) {
		    for (int j=0; j<seqB.length()+1; j++) {
			D[i][j] = (Mutation_FSM)fsmType.clone();
			Mutation_FSM.TraceBack_Info c = (Mutation_FSM.TraceBack_Info)D[i][j];
			c.set_tbdata(new Mutation_FSM.TraceBack_Data(i,j));
		    }
		}
		final_cell = (Mutation_FSM)fsmType.clone();
		((Mutation_FSM.TraceBack_Info)final_cell).set_tbdata(new Mutation_FSM.TraceBack_Data(-1,-1));
	    }

	    // Initialise the first cell of the DPA matrix
	    cell(D,0,0).init_counts(initialCounts); // Initialise counts
	    if (localAlign)
		cell(D,0,0).init_val(Double.POSITIVE_INFINITY); 
	    else
		cell(D,0,0).init_val(0); 

	    final_cell.init_val(Double.POSITIVE_INFINITY);

	    if (verbose>=1) { 
		System.out.println("\n\nIteration: " + iter);
		System.out.println(model);
		System.out.println(cell(D,0,0).paramsToString());
		//				   (localAlign ? "PalignChar="+PalignChar : ""));
	    }

	    // Do the DPA!
	    // NB. The Mutation_FSM calc function calculates the value of this cell
	    //     on the three neighbouring cell.  This is the reverse of the way the
	    //     DPA is usually expressed.
	    for (int i=0; i<seqA.length()+1; i++) {

		// Reset the next row of the DPA matrix
		if (!doTraceBack)
		    for (int j=0; j<seqB.length()+1; j++) {
			cell(D,i+1,j).reset();
		    }

		for (int j=0; j<seqB.length()+1; j++) {
		    Mutation_FSM v = (i==seqA.length() ? null : cell(D,i+1,j));
		    Mutation_FSM h = (j==seqB.length() ? null : cell(D,i,j+1));
		    Mutation_FSM d = (i==seqA.length() || j==seqB.length() ? null : cell(D,i+1,j+1));
		    char aChar = (i==seqA.length() ? '-' : seqA.charAt(i));
		    char bChar = (j==seqB.length() ? '-' : seqB.charAt(j));

		    if (localAlign) {
			double val;

			// Compute contribution of a local alignment that starts at (i,j)
			val = modelA.encodeCumulative(i) +  modelB.encodeCumulative(j);
			//val += encode_length(i);
			//val += encode_length(j);
			cell(D,i,j).or(val, initialCounts);


			// Compute contribution of a local alignment that ends at (i,j)
			val = cell(D,i,j).get_val() +
			    (modelA.encodeCumulative( seqA.length() ) - modelA.encodeCumulative(i)) +
			    (modelB.encodeCumulative( seqB.length() ) - modelB.encodeCumulative(j));

			//val += encode_length(seqA.length()-i);
			//val += encode_length(seqB.length()-j);

			//val += -MyMath.log2(1-PalignChar); // No more alignment characters
			//val += -MyMath.log2(1-PcharEncode); // No more alignment characters

			if (doTraceBack)
			    ((Mutation_FSM.TraceBack_Info)final_cell).or(val, 
									 cell(D,i,j).get_counts(),
									 cell(D,i,j));
			else
			    final_cell.or(val, cell(D,i,j).get_counts());


			//cell(D,i,j).add(-MyMath.log2(PalignChar), 0); // Count 0 is count for Palign
			//cell(D,i,j).add(-MyMath.log2(PcharEncode), 0); // Count 0 is count for Palign
		    }

		    cell(D,i,j).calc(h, v, d, aChar, bChar, i, j);

		    if (verbose>=3) {
			System.err.println("Calc outputs from D["+i+"]["+j+"]");
			System.err.println(cell(D,i,j));
		    }

		}
	    }


	    if (!localAlign)
		final_cell = cell(D, seqA.length(),seqB.length());
	    else {
		// Must be some _ends_ if there it is local alignment.
		/*
		if (!doTraceBack) 
		    final_cell.or(cell(D, seqA.length(), seqB.length()).get_val(), 
				  cell(D, seqA.length(), seqB.length()).get_counts());
		else
		    ((Mutation_FSM.TraceBack_Info)
		     final_cell).or(cell(D, seqA.length(), seqB.length()).get_val(), 
				    cell(D, seqA.length(), seqB.length()).get_counts(),
				    cell(D, seqA.length(), seqB.length()));
		*/
	    }
	    
	    double encAlignModel = final_cell.encode_params();
	    double encAlignment = encAlignModel + final_cell.get_val();
	    if (localAlign) {
		// Add a cost for the length of the alignment.
		// Assume we know the length of the sequences. Need to encode the start and end
		// of the alignment.  Assume uniform over all positions.  Have 4 cut-points to
		// encode, but only need to encode 3 (kind of).
		// So encode 2 from the shorter sequence, and 1 from the longer.
		double l1 = (seqA.length() < seqB.length() ? seqA.length() : seqB.length());
		double l2 = (seqA.length() > seqB.length() ? seqA.length() : seqB.length());
		encAlignment += MyMath.log2(l1) * 2 - 1;
		encAlignment += MyMath.log2(l2);
	    }


	    double encA = modelA.encodeCumulative( seqA.length() );
	    double encB = modelB.encodeCumulative( seqB.length() );
	    double encNull = encA + encB;

	    if (localAlign) {
		// Encode lengths for the null theory
		// Note that global alignments ignore lengths, ignore for null when doing global.
		//		encNull += encode_length(seqA.length());
		//		encNull += encode_length(seqB.length());
	    }

	    if (doTraceBack) {
		String s = getAlignment(D, (Mutation_FSM.TraceBack_Info)final_cell);
		System.out.println("ALIGNMENT:\n"+s);
	    }

	    // Get new parameters for next iteration
	    p = fsmType.counts_to_params(final_cell.get_counts());
	    /*
	    if (localAlign) {
		double n = final_cell.get_counts().get(0);
		p.put("PalignChar", (n-1)/n);
	    }
	    */

	    if (verbose>=2) {
		System.out.println("encA="+encA+" encB="+encB+" encNull="+(encNull));

		System.out.print("Mutual Encoding = " + encAlignment + " bits");
		System.out.println(" (model=" + encAlignModel + " data="+(encAlignment-encAlignModel)+")");
		System.out.println("\nCounts:\n"+final_cell.get_counts());
		System.out.println("\nEstimated Parameters:\n"+ p);
	    }
	    System.out.println((encAlignment < encNull ? 
				"related" : "unrelated") + " ("+(encNull-encAlignment)+")" +
			       "  log odds ratio = "+(encNull-encAlignment)+" bits");


	    if (bestDiff < encNull-encAlignment) bestDiff = encNull-encAlignment;

	    if (iter>0 && verbose>=1 && encAlignment>lastAlignment) {
		System.err.println("NON-CONVERGENCE: this="+encAlignment+" last="+lastAlignment);
	    }

	    // Done we little change in alignment length
	    if (iter>0 && lastAlignment-encAlignment < 0.05)
		break;

	    lastAlignment = encAlignment;
	    iter++;
	} // End iterations

	return bestDiff;
    }

    public static void main(String args[]) {
	CommandLine cmdLine = new CommandLine();
	cmdLine.addInt("markov", -1, "Order of Markov Model to use for sequence models.");
	cmdLine.addInt("iterations", -1, "Maximum number of iterations.");
	cmdLine.addBoolean("linearCosts", true, "Use linear gap costs.");
	cmdLine.addBoolean("sumAlignments", true, "Sum over all alignments.");
	cmdLine.addBoolean("local", false, "Compute using local alignments.");
	cmdLine.addInt("verbose", 0, "Display verbose output (larger num means more verbosity).");
	cmdLine.addString("params", "", "Params to pass to all classes (comma separated)");

	args = cmdLine.parseLine(args);

	AlignCompress a = new AlignCompress();

	if (args==null || args.length > 2) {
	    System.err.println("Usage: java AlignCompress [options]\n" + cmdLine.usage());
	    System.err.println("       sequences can be provided on the commandline, or from stdin");
	    System.exit(1);
	} else if (args.length == 2) {
	    // Two strings on commandline, assume they are sequences
	    a.seqA = args[0];
	    a.seqB = args[1];
	} else if (args.length == 1) {
	    // One arg, assume it is a filename to read the sequences from
	    try {
		BufferedReader in = new BufferedReader(new FileReader(args[0]));
		a.seqA = in.readLine();
		a.seqB = in.readLine();
	    } catch (Exception e) {
		System.err.println("Error reading '"+args[0]+"': "+e);
	    }
	} else if (args.length == 0) {
	    // No args, read sequence from stdin.
	    try {
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		a.seqA = in.readLine();
		a.seqB = in.readLine();
	    } catch (Exception e) {
		System.err.println("Error stdin: "+e);
	    }
	}

	if (a.seqA==null || a.seqB==null) {
	    System.err.println("Unable to read both sequences");
	    System.exit(1);
	}

	a.markovOrder     = cmdLine.getIntVal("markov");
	a.maxIterations   = cmdLine.getIntVal("iterations");
	a.linearCosts     = cmdLine.getBooleanVal("linearCosts");
	a.sumAlignments   = cmdLine.getBooleanVal("sumAlignments");
	a.localAlign      = cmdLine.getBooleanVal("local");
	a.verbose         = cmdLine.getIntVal("verbose");
	a.paramString     = cmdLine.getStringVal("params");

	a.alphabet = new char[] {'a', 't', 'g', 'c'};

	System.out.println("-log odds ratio = " + a.doAlign() + " bits");

    } // End main()
}
