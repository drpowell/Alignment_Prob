/*
 * Copyright (c) David Powell <david@drp.id.au>
 * 
 * 
 * This file is part of FuzzyLZ
 * 
 * FuzzyLZ is a program orginally intended for the compression of DNA sequeces.
 * It can be viewed as a compression model like Lempel-Ziv 77, but instead of
 * exact matches, allowing matches that contain inserts/deletes/mismatches.
 *  
 */

package fuzzyLZ;

import java.io.*;

import common.*;

class Plot implements Serializable {
    byte red[][], grn[][], blu[][];

    double scale;

    int numRows, numCols;

    int startRow;

    Plot(int rows, int columns, int maxColumns, int maxRows) {
        this(rows, columns, maxColumns, maxRows, 0);
    }

    Plot(int rows, int columns, int maxColumns, int maxRows, int startRow) {
        this.startRow = startRow;

        scale = MyMath.min2(1.0 * MyMath.min2(columns, maxColumns) / columns,
                1.0 * MyMath.min2(rows - startRow, maxRows) / rows);

        numRows = (int) (scale * (rows - startRow));
        numCols = (int) (scale * columns);
        red = new byte[numRows + 1][numCols+1];
        grn = new byte[numRows + 1][numCols+1];
        blu = new byte[numRows + 1][numCols+1];

    }

    void put(int row, int col, double r, double g, double b) {
        row = (int) (scale * (row - startRow));
        col = (int) (scale * col);
        red[row][col] = (byte) (r * 255);
        grn[row][col] = (byte) (g * 255);
        blu[row][col] = (byte) (b * 255);
    }

    void putMax(int row, int col, double r, double g, double b) {
        row = (int) (scale * (row - startRow));
        col = (int) (scale * col);

        byte r1 = (byte) MyMath.max2(red[row][col], r * 255);
        byte g1 = (byte) MyMath.max2(grn[row][col], g * 255);
        byte b1 = (byte) MyMath.max2(blu[row][col], b * 255);

        red[row][col] = r1;
        grn[row][col] = g1;
        blu[row][col] = b1;
    }

    void save(String fname, String comments) {
        if (FuzzyLZ.DEBUG >= 1)
            System.out.println("Writing image '" + fname + "'");
        try {
            File f = new File(fname);
            DataOutputStream out = new DataOutputStream(
                    new BufferedOutputStream(new FileOutputStream(f)));

            out.writeBytes("P6\n");
            out.writeBytes("#Created by FuzzyLZ\n");
            if (comments.length() > 0)
                out.writeBytes("#" + comments + "\n");
            out.writeBytes(numCols + " " + numRows + "\n");
            out.writeBytes("255\n");

            for (int r = 0; r < numRows; r++) {
                for (int c = 0; c < numCols; c++) {
                    out.writeByte(red[r][c]);
                    out.writeByte(grn[r][c]);
                    out.writeByte(blu[r][c]);
                }
            }

            out.close();
        } catch (Exception e) {
            System.err.println("Error writing file: " + e);
        }
        if (FuzzyLZ.DEBUG >= 1)
            System.out.println("Done image");
    }

    public static void main(String args[]) {
        Plot p = new Plot(100, 100, 50, 50);
        for (int i = 0; i < 100; i++) {
            p.put(i, i, 0, 1, 0);
            p.put(i, 99 - i, 0, 0, 1);
            p.put(50, i, 1, 0, 0);
            p.put(i, 50, 1, 0, 0);
        }
        p.save("out.ppm", "");
    }
}
