
package alignCompress;

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

    int markovOrder;
    int numIterations;
    int verbose;
    boolean linearCosts;
    boolean localAlign;

    public AlignCompress() {
    }

    public double doAlign() {
	System.err.println("# SeqA = "+seqA+"\n# SeqB = "+seqB+
			   "\n# markov order="+markovOrder+
			   "\n# num iterations="+numIterations+
			   "\n# linearCosts="+linearCosts+
			   "\n# local Alignment="+localAlign+
			   "\n# verbosity="+verbose+
			   "\n");

	Params p = new Params();


	BufferModel modelA = new BufferModel( new MarkovN(markovOrder, alphabet), seqA, alphabet );
	BufferModel modelB = new BufferModel( new MarkovN(markovOrder, alphabet), seqB, alphabet );

	//System.err.println("modelA\n"+modelA+"\n");
	//System.err.println("modelB\n"+modelB+"\n");

	int mdlCounts = Model_SeqAB.required_counts();
	int fsmCounts = (linearCosts ? FSM_3State.required_counts() : FSM_Prob.required_counts()); 
	int totCounts = mdlCounts + fsmCounts;

	double bestDiff = -Double.POSITIVE_INFINITY;

	for (int iter=0; iter<numIterations; iter++) {
	    // EncNonAlignChar - (local alignments) encode that a non-aligned character follows
	    // EncEndNonAlign  - (local alignments) next will be the alignment
	    // EncAlignChar - (local alignments) encode that an alignment pair follows
	    // EncEndAlign  - (local alignments) encode no more of the alignment to follow
	    //        The use of these corresponds to encoding the lengths using a geometric distribution.
	    //        This 'aint ideal but it's a start.  The values here seem reasonable.
	    //                                         (should they be estimated as well?)
	    double EncNonAlignChar, EncEndNonAlign;
	    double EncAlignChar, EncEndAlign;

	    {
		EncNonAlignChar = -(MyMath.log2(seqA.length()+seqB.length()) - 
				    MyMath.log2(seqA.length()+seqB.length() + 1));
		EncEndNonAlign  = -(MyMath.log2(1) -
				    MyMath.log2(seqA.length()+seqB.length() + 1));

		EncAlignChar = -(MyMath.log2(seqA.length()+seqB.length()) - 
				 MyMath.log2(seqA.length()+seqB.length() + 2));
		EncEndAlign  = -(MyMath.log2(2) -
				 MyMath.log2(seqA.length()+seqB.length() + 2));
	    }

	    int countPos = 0;
	    Two_Seq_Model_Counts model = new Model_SeqAB(p, modelA, modelB, countPos);
	    countPos += mdlCounts;
	    Mutation_FSM.With_Counts fsmType;
	    if (linearCosts)
		fsmType = new FSM_3State(model, p, totCounts, countPos);
	    else
		fsmType = new FSM_Prob(model, p, totCounts, countPos);
	    countPos += fsmCounts;
	    Misc.assert(countPos == totCounts, "Internal error: countPos!=totCounts");

	    Counts emptyCounts = new Counts(totCounts);
	    for (int i=0; i<totCounts; i++) 
		emptyCounts.inc(i, 1); // Start all counts at 1.
	    
	    
	    // Setup the DPA matrix
	    Mutation_FSM.With_Counts D[][] = new Mutation_FSM.With_Counts[2][seqB.length()+1];
	    for (int i=0; i<2; i++) {
		for (int j=0; j<seqB.length()+1; j++) {
		    D[i][j] = (Mutation_FSM.With_Counts)fsmType.clone();
		}
	    }

	    // Initialise the first cell of the DPA matrix
	    D[0][0].init_val(0); 
	    D[0][0].init_counts(emptyCounts); // Initialise counts

	    if (verbose>=1) { 
		System.out.println("Iteration: " + iter);
		System.out.println(model);
		System.out.println(D[0][0].paramsToString()); 
	    }

	    // Do the DPA!
	    // NB. The Mutation_FSM calc function calculates the value of this cell
	    //     on the three neighbouring cell.  This is the reverse of the way the
	    //     DPA is usually expressed.
	    for (int i=0; i<seqA.length()+1; i++) {

		// Reset the next row of the DPA matrix
		for (int j=0; j<seqB.length()+1; j++) {
		    D[(i+1)%2][j].reset();
		}

		for (int j=0; j<seqB.length()+1; j++) {
		    Mutation_FSM v = (i==seqA.length() ? null : D[(i+1)%2][j]);
		    Mutation_FSM h = (j==seqB.length() ? null : D[ i   %2][j+1]);
		    Mutation_FSM d = (i==seqA.length() || j==seqB.length() ? null : D[(i+1)%2][j+1]);
		    char aChar = (i==seqA.length() ? '-' : seqA.charAt(i));
		    char bChar = (j==seqB.length() ? '-' : seqB.charAt(j));

		    if (localAlign) {
			// Compute contribution of a local alignment that starts at (i,j)
			double val = modelA.encodeCumulative(i) +  modelB.encodeCumulative(j);
			val += i*EncNonAlignChar + EncEndNonAlign;
			val += j*EncNonAlignChar + EncEndNonAlign;
			D[i%2][j].or(val, emptyCounts);

			// Encode that an alignment character is coming next
			D[i%2][j].add(EncAlignChar, -1);
		    }

		    if (verbose>=3) {
			System.err.println("Calc outputs from D["+i+"]["+j+"]");
			System.err.println(D[i%2][j]);
		    }

		    
		    D[i%2][j].calc(h, v, d, aChar, bChar, i, j);

		    if (localAlign) {
			// Compute contribution of a local alignment that ends at (i,j)
			double val = D[i%2][j].get_val() +
			    (modelA.encodeCumulative( seqA.length() ) - modelA.encodeCumulative(i)) +
			    (modelB.encodeCumulative( seqB.length() ) - modelB.encodeCumulative(j));

			// Remove encoding that an alignment char is next.
			// Then encode that there is no more alignment chars.
			val -= EncAlignChar;
			val += EncEndAlign;

			val += (seqA.length()-i)*EncNonAlignChar + EncEndNonAlign;
			val += (seqB.length()-j)*EncNonAlignChar + EncEndNonAlign;
			D[seqA.length()%2][seqB.length()].or(val, D[i%2][j].get_counts());
		    }
		}
	    }

	    double encAlignModel = D[seqA.length()%2][seqB.length()].encode_params();
	    double encAlignment = encAlignModel + D[seqA.length()%2][seqB.length()].get_val();

	    double encA = modelA.encodeCumulative( seqA.length() );
	    double encB = modelB.encodeCumulative( seqB.length() );
	    double encNull = encA + encB;

	    if (localAlign) encNull += (seqA.length()+seqB.length())*EncNonAlignChar + EncEndNonAlign;

	    if (bestDiff < encNull-encAlignment) bestDiff = encNull-encAlignment;

	    // Get new parameters for next iteration
	    p = fsmType.counts_to_params(D[seqA.length()%2][seqB.length()].get_counts());

	    if (verbose>=2) {
		System.out.println("encA="+encA+" encB="+encB+" encNull="+(encNull));

		System.out.print("Mutual Encoding = " + encAlignment + " bits");
		System.out.println(" (model=" + encAlignModel + " data="+(encAlignment-encAlignModel)+")");
		System.out.println((encAlignment < encNull ? 
				    "related" : "unrelated") + " ("+(encNull-encAlignment)+")");
	    
		System.out.println("\nCounts:\n"+D[seqA.length()%2][seqB.length()].get_counts());
		System.out.println("\nEstimated Parameters:\n"+
				   fsmType.counts_to_params(D[seqA.length()%2][seqB.length()].get_counts()));
	    }
	} // End iterations

	return bestDiff;
    }

    public static void main(String args[]) {
	CommandLine cmdLine = new CommandLine();
	cmdLine.addInt("markov", -1, "Order of Markov Model to use for sequence models.");
	cmdLine.addInt("iterations", 1, "Number of iterations.");
	cmdLine.addBoolean("linearCosts", true, "Use linear gap costs.");
	cmdLine.addBoolean("local", false, "Compute using local alignments.");
	cmdLine.addInt("verbose", 0, "Display verbose output (larger num means more verbosity).");

	args = cmdLine.parseLine(args);
	if (args == null || args.length != 2) {
	    System.err.println("Usage: java AlignCompress [options] seqA seqB\n" + cmdLine.usage());
	    System.exit(1);
	}

	AlignCompress a = new AlignCompress();

	a.markovOrder     = cmdLine.getIntVal("markov");
	a.numIterations   = cmdLine.getIntVal("iterations");
	a.linearCosts = cmdLine.getBooleanVal("linearCosts");
	a.localAlign  = cmdLine.getBooleanVal("local");
	a.verbose     = cmdLine.getIntVal("verbose");

	a.seqA = args[0];
	a.seqB = args[1];
	a.alphabet = new char[] {'a', 't', 'g', 'c'};

	System.out.println("-log odds ratio = " + a.doAlign() + " bits");

    } // End main()
}
