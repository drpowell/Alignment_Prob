
package common;

public final class Params {
    private String names[];
    private double vals[];
    private int num,room;

    private void alloc(int r) {
	if (room >= r)
	    return;

	String[] n = new String[r];
	double[]  v = new double[r];

	if (room>0) {		// Copy over any data
	    for (int i=0; i<num; i++) {
		n[i] = names[i];
		v[i] = vals[i];
	    }
	}

	room = r;
	names = n;
	vals = v;
    }

    public Params() {
	num=0;
	room=0;
	alloc(10);
    }

    public Params put(String s, double v) {
	int id = get_id(s);
        if (id<0) {
	    if (num==room) alloc(room*2);
	    names[num] = new String(s);
	    id = num;
	    num++;
	}
        vals[id] = v;

	return this;
    }

    private int get_id(String s) {
	for (int i=0; i<num; i++) {
	    if (s.equalsIgnoreCase(names[i]))
		return i;
	}
	return -1;
    }

    public boolean exists(String s) {
	return get_id(s)>=0;
    }
    
    public double get(String s) {
	int id = get_id(s);
	Misc.assert(id>=0, "Attempt Params.get() with non-existent key");
	return vals[id];
    }

    public int get_num() { return num; }

    public String get_name_by_id(int id) { return names[id]; }

    
    public void join(Params p) {
	for (int i=0; i<p.num; i++) {
	    put(p.names[i], p.vals[i]);
	}
    }

    public String toString() {
	StringBuffer r = new StringBuffer();
	for (int i=0; i<num; i++) {
	    r.append(names[i]+"="+vals[i]+"\n");
	}
	return r.toString();
    }

}

