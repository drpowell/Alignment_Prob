
// DNA - read DNA sequence from a file in either raw format
// or genbank format.  Note, it guesses the format in a
// pretty stupid manner.

package common;

import java.io.*;
import java.util.*;
import java.util.zip.*;

public class DNA {
    BufferedReader in;
    String filename;
    char[] seq;

    DNA(String fname) {
	filename = fname;
	in = null;
    }

    DNA(BufferedReader r) {
	in = r;
	filename = null;
    }

    public char[] sequence() {
	return seq;
    }

    public String toString() {
	return new String(seq);
    }

    protected static BufferedReader openFile(String fname) {
	BufferedReader in = null;
	try {
	    int magic;
	    InputStream is = new BufferedInputStream(new FileInputStream(fname));

	    is.mark(4);
	    magic = is.read();
	    magic += is.read()*256;
	    is.reset();

	    if (magic == GZIPInputStream.GZIP_MAGIC)
		is = new GZIPInputStream(is);

	    in = new BufferedReader(new InputStreamReader(is));
	} catch (IOException e) {
	    System.err.println("Error reading '"+fname+"' "+e);
	}
	return in;
    }

    public static DNA guess_format(String filename) {
	BufferedReader in = openFile(filename);
	if (in==null) return null;

	try {
	    in.mark(10);
	    char[] buf = new char[10];
	    in.read(buf, 0, 10);
	    in.reset();

	    if (new String(buf).startsWith("LOCUS"))
		return new DNA.Genbank(in);
	    else if (new String(buf).startsWith(">"))
		return new DNA.FASTA(in);
	    else
		return new DNA.Raw(in);

	} catch (IOException e) {
	    System.err.println("Error reading '"+filename+"' "+e);
	}
	return null;
    }




    static class Raw extends DNA {
	Raw(BufferedReader r) {
	    super(r);
	    read(r);
	}

	Raw(String fname) {
	    super(fname);
	    
	    read(openFile(fname));
	}

	private void read(BufferedReader in) {
	    System.err.println("Try to read as raw");
	    try {
		StringBuffer seqBuf = new StringBuffer();
		char[] buf = new char[1024];
		int n;
		while( (n=in.read(buf))>=0) {
		    for (int i=0; i<n; i++) {
			char c = buf[i];
			if (Character.isWhitespace(c) || Character.isDigit(c)) continue;	
			seqBuf.append(c);
		    }
		}
		seq = new char[seqBuf.length()];
		seqBuf.getChars(0, seqBuf.length(), seq, 0);
	    } catch (IOException e) {
		System.err.println("Error reading '"+filename+"' "+e);
	    }
	}
    }

    static class FASTA extends DNA {
	String name;
	FASTA(BufferedReader r) {
	    super(r);
	    name = null;

	    read(r);
	}

	FASTA(String fname) {
	    super(fname);
	    name = null;
	    
	    read(openFile(fname));
	}

	private void read(BufferedReader in) {
	    System.err.println("Try to read as FASTA");
	    try {
		StringBuffer seqBuf   = new StringBuffer();
		String line;
		while ( (line=in.readLine()) != null) {
		    if (name==null && line.startsWith(">")) {
			name = line.substring(1);
			continue;
		    }

		    if (name==null) throw(new IOException("Bad FASTA format. Name not found"));

		    for (int i=0; i<line.length(); i++) {
			char c = line.charAt(i);
			
			if (Character.isWhitespace(c) || Character.isDigit(c)) continue;
			seqBuf.append(c);
		    }
		}
		seq = new char[seqBuf.length()];
		seqBuf.getChars(0, seqBuf.length(), seq, 0);
	    } catch (IOException e) {
		System.err.println("Error reading '"+filename+"' "+e);
	    }	
	}
    }

    static class Genbank extends DNA {
	String name;
	StringBuffer other;
	StringBuffer features;

	Genbank(BufferedReader r) {
	    super(r);
	    name     = null;
	    other    = null;
	    features = null;

	    read(r);
	}

	Genbank(String fname) {
	    super(fname);
	    name     = null;
	    other    = null;
	    features = null;

	    read(openFile(fname));
	}

	private void read(BufferedReader in) {
	    System.err.println("Try to read as genbank");
	    try {
		other    = new StringBuffer();
		features = new StringBuffer();
		StringBuffer seqBuf   = new StringBuffer();
		int state = 0;
		String line;
		while ( (line=in.readLine()) != null) {
		    if (line.startsWith("LOCUS")) {
			name = line.substring(5);
			state = 1;
			continue;
		    }

		    if (state==0 || line.length()==0) continue;

		    if (line.startsWith("ORIGIN")) {
			state = 2;
			continue;
		    }
		    
		    if (line.startsWith("FEATURES")) {
			state = 3;
			continue;
		    }

		    if (state == 2) {
			if (line.startsWith("//")) {
			    state = 1;
			    continue;
			}

			for (int i=0; i<line.length(); i++) {
			    char c = line.charAt(i);
			    if (Character.isWhitespace(c) || Character.isDigit(c)) continue;
			    seqBuf.append(c);
			}
			continue;
		    }

		    if (!Character.isWhitespace(line.charAt(0)))
			state = 1;

		    if (state == 3) {
			features.append(line + "\n");
			continue;
		    }

		    other.append(line + "\n");
		}

		seq = new char[seqBuf.length()];
		seqBuf.getChars(0, seqBuf.length(), seq, 0);
	    } catch (IOException e) {
		System.err.println("Error reading '"+filename+"' "+e);
	    }
	}

    }

    public static void main(String[] args) {
	DNA d = DNA.guess_format(args[0]);
	System.out.println(d.sequence());
    }			       
}

