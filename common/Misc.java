
package common;

import java.io.*;

public final class Misc {
    public static void assert(boolean e, String s) {
        if (!e) error(s);
    }

    public static void error(String s) {
        System.err.println("ERROR: " + s);

	Misc o = (Misc)(new Object()); // Force a crash.  want stack trace.
        System.exit(1);
    }

    static String readString(InputStream in) {
	StringBuffer s = new StringBuffer();
        try {
            byte[] buf = new byte[1024];
            int n;
            while( (n=in.read(buf))>=0) {
                s.append( new String(buf,0,n) );
            }
        } catch (Exception e) {
            System.err.println("Unable to read stdin");
            System.exit(1);
        }
	return s.toString();
    }
}

