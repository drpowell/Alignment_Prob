
package common;

public class CommandLine {
    int maxOptions = 100;
    String[]  opts;
    char[]    types;
    Object[]  defaults;
    String[]  helps;
    boolean[] option_set;
    int numOptions;

    Object[] values;

    public CommandLine() {
	opts     = new String[maxOptions];
	types    = new char  [maxOptions];
	defaults = new Object[maxOptions];
	helps    = new String[maxOptions];

	option_set = new boolean[maxOptions];
	values   = new Object[maxOptions];
	numOptions = 0;
    }

    public String usage() {
	StringBuffer res = new StringBuffer();
	for (int i=0; i<numOptions; i++) {
	    res.append("  --" + opts[i]);

	    int len = opts[i].length();
	    switch (types[i]) {
	    case 'b': break;
	    case 'i':
		res.append("=i");
		len+=2;
		break;
	    case 'f':
		res.append("=d");
		len+=2;
		break;
	    case 's':
		res.append("=s");
		len+=2;
		break;
	    }

	    for (int j=0; j<20-len; j++)
		res.append(" ");
	    res.append(helps[i]);
	    //res.append(" (default='"+defaults[i]+"')");
	    res.append("\n                        (default='"+defaults[i]+"')");
	    res.append("\n");
	}

	return res.toString();
    }

    private void addOption(String opt, char type, Object def, String help) {
	if (numOptions >= maxOptions) {
	    System.err.println("Too many options to CommandLine class");
	    System.exit(1);
	}
	opts[numOptions]       = opt;
	types[numOptions]      = type;
	defaults[numOptions]   = def;
	helps[numOptions]      = help; 
	option_set[numOptions] = false;
	numOptions++;
    }

    public void addBoolean(String opt, boolean def, String help) {	
	addOption(opt, 'b', new Boolean(def), help); 
    }

    public void addInt(String opt, int def, String help) { 
	addOption(opt, 'i', new Integer(def), help); 
    }

    public void addDouble(String opt, double def, String help) { 
	addOption(opt, 'f', new Double(def), help); 
    }

    public void addString(String opt, String def, String help) { 
	addOption(opt, 's', def, help); 
    }

    public boolean optionSet(String opt) {
	int o = isOption(opt, 1);
	if (o<0) {
	    System.err.println("Attempt to lookup non-defined option '"+opt+"'");
	    return false;
	}
	return option_set[o];
    }

    Object getVal(String opt) {
	int o = isOption(opt, 1);
	if (o<0) {
	    System.err.println("Attempt to lookup non-defined option '"+opt+"'");
	    return null;
	}
	return values[o];
    }

    public int getIntVal(String opt) {
	int o = isOption(opt, 1);
	if (o<0) {
	    System.err.println("Attempt to lookup non-defined option '"+opt+"'");
	    return 0;
	}

	if (types[o] != 'i') {
	    System.err.println("Option '"+opt+"' is not an int option in getIntVal");
	    return 0;
	}
	
	return ((Integer)values[o]).intValue();
    }

    public double getDoubleVal(String opt) {
	int o = isOption(opt, 1);
	if (o<0) {
	    System.err.println("Attempt to lookup non-defined option '"+opt+"'");
	    return 0;
	}

	if (types[o] != 'f') {
	    System.err.println("Option '"+opt+"' is not a double option in getDoubleVal");
	    return 0;
	}
	
	return ((Double)values[o]).doubleValue();
    }

    public String getStringVal(String opt) {
	int o = isOption(opt, 1);
	if (o<0) {
	    System.err.println("Attempt to lookup non-defined option '"+opt+"'");
	    return null;
	}

	if (types[o] != 's') {
	    System.err.println("Option '"+opt+"' is not a string option in getStringVal");
	    return null;
	}
	
	return (String)values[o];
    }

    public boolean getBooleanVal(String opt) {
	int o = isOption(opt, 1);
	if (o<0) {
	    System.err.println("Attempt to lookup non-defined option '"+opt+"'");
	    return false;
	}

	if (types[o] != 'b') {
	    System.err.println("Option '"+opt+"' is not a boolean option in getBooleanVal");
	    return false;
	}
	
	return ((Boolean)values[o]).booleanValue();
    }

    int isOption(String opt, int noDashOk) {
	if      (opt.startsWith("--")) opt = opt.substring(2);
	else if (opt.startsWith("-"))  opt = opt.substring(1);
	else if (noDashOk == 0)
	    return -1;

	if (opt.indexOf("=")>=0) {
	    opt = opt.substring(0, opt.indexOf("="));
	}

	int match = -2;
	for (int i=0; i<numOptions; i++) {
	    if (opt.length() > opts[i].length())
		continue;

	    String optStr = opts[i].substring(0, opt.length());
	    if ( opt.compareToIgnoreCase(optStr)==0 ) {
		if (match>=0) {
		    System.err.println("Ambiguous option '"+opt+"' could be '"+
				       opts[i]+"' or '"+opts[match]+"'");
		    return match;
		}
		match = i;
	    }
	}
	return match;
    }

    public String[] parseLine(String[] args) {
	int i=0;
	int j=0 ;
	String[] res = new String[args.length];

	// First setup defaults
	for (int k=0; k<numOptions; k++)
	    values[k] = defaults[k];

	while (i<args.length) {
	    int o = isOption(args[i], 0);
	    if (o==-2) {
		//System.err.println("Unknown option '"+args[i]+"'");
		return null;
	    }

	    if (o>=0) {
		option_set[o] = true;

		try {
		    switch (types[o]) {
		    case 'b':
			if (args[i].indexOf("=") >=0 ) {
			    String s = args[i].substring(args[i].indexOf("=")+1);
			    if (s.equalsIgnoreCase("true"))
				values[o] = Boolean.TRUE;
			    else if (s.equalsIgnoreCase("yes"))
				values[o] = Boolean.TRUE;
			    else if (s.equalsIgnoreCase("on"))
				values[o] = Boolean.TRUE;
			    else if (s.equalsIgnoreCase("1"))
				values[o] = Boolean.TRUE;
			    else if (s.equalsIgnoreCase("false"))
				values[o] = Boolean.FALSE;
			    else if (s.equalsIgnoreCase("no"))
				values[o] = Boolean.FALSE;
			    else if (s.equalsIgnoreCase("off"))
				values[o] = Boolean.FALSE;
			    else if (s.equalsIgnoreCase("0"))
				values[o] = Boolean.FALSE;
			    else {
				System.err.println("Unknown boolean option parameter '"+s+"'");
				return null;
			    }
			} else
			    values[o] = Boolean.TRUE;
			break;
		    case 'i':
			if (args[i].indexOf("=") >=0 ) 
			    values[o] = new Integer(args[i].substring(args[i].indexOf("=")+1));
			else {
			    values[o] = new Integer(args[i+1]);
			    i++;
			}
			break;
		    case 'f':
			if (args[i].indexOf("=") >=0 ) 
			    values[o] = new Double(args[i].substring(args[i].indexOf("=")+1));
			else {
			    values[o] = new Double(args[i+1]);
			    i++;
			}
			break;
		    case 's':
			if (args[i].indexOf("=") >=0 ) 
			    values[o] = args[i].substring(args[i].indexOf("=")+1);
			else {
			    values[o] = args[i+1];
			    i++;
			}
			break;
		    }
		} catch (Exception e) {
		    System.err.println("Bad option '"+args[i]+"'");
		    return null;
		}
	    } else {
		res[j++] = args[i];
	    }

	    i++;
	}

	String[] r = new String[j];
	for (int k=0; k<j; k++) r[k]=res[k];
	return r;
    }

    public static void main(String args[]) {
	CommandLine cmdLine = new CommandLine();
	cmdLine.addInt("foo", 7, "Some integer");
	cmdLine.addBoolean("fom", false, "A bool");
	cmdLine.addDouble("doub", 8.42, "A double");
	cmdLine.addString("str", "cthulu", "A string");

	args = cmdLine.parseLine(args);
	if (args == null) {
	    System.err.println("Usage: java CommandLine [options]\n" + cmdLine.usage());
	    System.exit(1);
	}

	System.out.println("foo  = " + cmdLine.getVal("foo"));
	System.out.println("fom  = " + cmdLine.getVal("fom"));
	System.out.println("doub = " + cmdLine.getVal("doub"));
	System.out.println("str  = " + cmdLine.getVal("str"));
	System.out.println("Rest:");
	for (int i=0; i<args.length; i++) { System.out.println("args["+i+"] = "+args[i]); }
    }
}
