
package fuzzyLZ;

import java.io.*;

import common.*;

class Plot implements Serializable {
    int img[][];
    double xscale, yscale;
    int w,h;

    Plot(int width, int height, int maxWidth, int maxHeight) {
	w = MyMath.min2(width, maxWidth);
	h = MyMath.min2(height, maxHeight);
	img = new int[w][h];
	xscale = 1.0*w/width;
	yscale = 1.0*h/height;
    }

    void put(int x, int y, double r, double g, double b) {
	x = (int)(xscale*x);
	y = (int)(yscale*y);
	img[x][y] = ((byte)(r*255))<<16 | ((byte)(g*255))<<8 | ((byte)(b*255));
    }

    void putMax(int x, int y, double r, double g, double b) {
	x = (int)(xscale*x);
	y = (int)(yscale*y);

	byte r1 = (byte)MyMath.max2((img[x][y]>>16)&255, r*255);
	byte g1 = (byte)MyMath.max2((img[x][y]>> 8)&255, g*255);
	byte b1 = (byte)MyMath.max2((img[x][y]    )&255, b*255);
	
	img[x][y] = (r1)<<16 | (g1)<<8 | (b1);
    }

    void save(String fname, String comments) {
	System.err.println("Writing image");
	try {
	    File f = new File(fname);
	    DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f)));

	    out.writeBytes("P6\n");
	    out.writeBytes("#Created by FuzzyLZ\n");
	    if (comments.length()>0) 
		out.writeBytes("#"+comments+"\n");
	    out.writeBytes(w + " " + h + "\n");
	    out.writeBytes("255\n");

	    for (int i=0; i<w; i++) {
		for (int j=0; j<h; j++) {
		    out.writeByte( (img[i][j]>>16)&255 );
		    out.writeByte( (img[i][j]>> 8)&255 );
		    out.writeByte( (img[i][j]    )&255 );
		}
	    }

	    out.close();
	} catch (Exception e) {
	    System.err.println("Error writing file: "+e);
	}
	System.err.println("Done image");
    }



    public static void main(String args[]) {
	Plot p = new Plot(100,100,50,50);
	for (int i=0; i<100; i++) {
	    p.put(i,i, 0, 255, 0);
	    p.put(i,99-i, 0, 0, 255);
	    p.put(50, i, 255, 0, 0);
	    p.put(i, 50, 255, 0, 0);
	}
	p.save("out.ppm", "");
    }
}
